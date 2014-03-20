import re
import argparse
import cPickle as pickle

import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.image as afwImage
import lsst.coadd.utils as coaddUtils
from lsst.meas.algorithms import CoaddPsf
from lsst.pex.config import Config, Field, ConfigurableField, ConfigField, ListField
from lsst.pex.exceptions import LsstCppException, InvalidParameterException
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, SelectDataIdContainer, CoaddTaskRunner
from lsst.pipe.tasks.selectImages import BaseSelectImagesTask, WcsSelectImagesTask, BaseExposureInfo
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.tasks.assembleCoadd import (AssembleCoaddTask, AssembleCoaddConfig, _subBBoxIter,
                                           AssembleCoaddDataIdContainer,)
from hsc.pipe.tasks.processCoadd import SubaruProcessCoaddTask
from lsst.pipe.tasks.coaddHelpers import groupPatchExposures, getGroupDataRef
from lsst.pipe.base import Struct, DataIdContainer, ArgumentParser
from hsc.pipe.base.pbs import PbsPoolTask
from hsc.pipe.base.pool import Pool, abortOnError, NODE, Debugger
from hsc.pipe.tasks.background import BackgroundReferenceTask, MatchBackgroundsTask

from hsc.pipe.tasks.backgroundModels import Background, BackgroundConfig


#Debugger().enabled = True

DEBUGGING = False

class NullSelectImagesTask(BaseSelectImagesTask):
    """Select images by taking everything we're given without further examination

    This is useful if the examination (e.g., Wcs checking) has been performed
    previously, and we've been provided a good list.
    """
    def runDataRef(self, patchRef, coordList, makeDataRefList=True, selectDataList=[]):
        return Struct(
            dataRefList = [s.dataRef for s in selectDataList],
            exposureInfoList = [BaseExposureInfo(s.dataRef.dataId, None) for s in selectDataList],
            )


class SimpleAssembleCoaddConfig(AssembleCoaddConfig):
    matchBackgrounds = ConfigurableField(target=MatchBackgroundsTask, doc="Background matching")
    doBackgroundSubtraction = Field(dtype=bool, default=True, doc="Subtract background?")
    background = ConfigField(dtype=BackgroundConfig, doc="Background matching config")

    def setDefaults(self):
        super(SimpleAssembleCoaddConfig, self).setDefaults()
        self.badMaskPlanes += ["BAD", "CR", "EDGE", ]

class SimpleAssembleCoaddTask(AssembleCoaddTask):
    """Assemble a coadd from a set of coaddTempExp
    """
    ConfigClass = SimpleAssembleCoaddConfig

    def __init__(self, *args, **kwargs):
        super(SimpleAssembleCoaddTask, self).__init__(*args, **kwargs)
        self.makeSubtask("interpImage")
        self.makeSubtask("matchBackgrounds")
        self.makeSubtask("scaleZeroPoint")
        self.debug = False

    def run(self, dataRef, selectDataList=[]):
        """Assemble a coadd from a set of coaddTempExp

        The coadd is computed as a mean with optional outlier rejection.

        assembleCoaddTask only works on the dataset type 'coaddTempExp', which are 'coadd temp exposures'.
        Each coaddTempExp is the size of a patch and contains data for one run, visit or
        (for a non-mosaic camera it will contain data for a single exposure).

        coaddTempExps, by default, will have backgrounds in them and will require
        config.doMatchBackgrounds = True. However, makeCoaddTempExp.py can optionally create background-
        subtracted coaddTempExps which can be coadded here by setting
        config.doMatchBackgrounds = False.

        @param dataRef: data reference for a coadd patch (of dataType 'Coadd')
        Used to access the following data products (depending on the config):
        - [in] self.config.coaddName + "Coadd_tempExp"
        - [in] self.config.coaddName + "Coadd_bgRef"
        - [out] self.config.coaddName + "Coadd"

        @return: a Struct with fields:
        - coaddExposure: coadd exposure
        """
        skyInfo = self.getSkyInfo(dataRef)
        calExpRefList = self.selectExposures(dataRef, skyInfo, selectDataList=selectDataList)
        if len(calExpRefList) == 0:
            self.log.warn("No exposures to coadd")
            return
        self.log.info("Coadding %d exposures" % len(calExpRefList))

        tempExpRefList = self.getTempExpRefList(dataRef, calExpRefList)
        inputData = self.prepareInputs(dataRef, tempExpRefList)
        tempExpRefList = inputData.tempExpRefList
        self.log.info("Found %d %s" % (len(inputData.tempExpRefList), self.getTempExpDatasetName()))
        if len(inputData.tempExpRefList) == 0:
            self.log.warn("No coadd temporary exposures found")
            return

        coaddExp = self.assemble(skyInfo, inputData.tempExpRefList, inputData.imageScalerList,
                                 inputData.weightList, inputData.bgModelList)

        if self.config.doInterp:
            fwhmPixels = self.config.interpFwhm / skyInfo.wcs.pixelScale().asArcseconds()
            self.interpImage.interpolateOnePlane(maskedImage=coaddExp.getMaskedImage(), planeName="EDGE",
                                                 fwhmPixels=fwhmPixels)

        if self.config.doBackgroundSubtraction:
            # XXX Should be a two-step process, with a detection stage in-between
            bgModel = self.subtractBackground(coaddExp)
        else:
            bgModel = None

        if self.config.doWrite:
            self.writeCoaddOutput(dataRef, coaddExp)
            if bgModel:
                bgImage = bgModel.getImageF(coaddExp.getMaskedImage().getBBox(afwImage.PARENT))
                self.writeCoaddOutput(dataRef, bgImage, "bg")

        return Struct(coaddExposure=coaddExp, background=bgModel)

    def getTempExpRefList(self, patchRef, calExpRefList):
        """Generate list of warp data references

        @param dataRef: Data reference for patch
        @param calExpRefList: List of data references for input calexps
        @return List of warp data references
        """
        butler = patchRef.getButler()
        groupData = groupPatchExposures(patchRef, calExpRefList, self.getCoaddDatasetName(),
                                        self.getTempExpDatasetName())
        tempExpRefList = [getGroupDataRef(butler, self.getTempExpDatasetName(), g, groupData.keys) for
                          g in groupData.groups.keys()]
        return tempExpRefList

    def prepareInputs(self, patchRef, refList, tract):
        """Prepare the input warps for coaddition

        This involves measuring weights and constructing image scalers
        for each of the inputs.

        @param refList: List of data references to tempExp
        @param tract: Tract on which we're coadding
        @return Struct:
        - tempExprefList: List of data references to tempExp
        - weightList: List of weightings
        - imageScalerList: List of image scalers
        """
        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(self.getBadPixelMask())
        statsCtrl.setNanSafe(True)

        # compute tempExpRefList: a list of tempExpRef that actually exist
        # and weightList: a list of the weight of the associated coadd tempExp
        # and imageScalerList: a list of scale factors for the associated coadd tempExp
        tempExpRefList = []
        weightList = []
        imageScalerList = []
        bgModelList = []
        if self.config.doMatchBackgrounds:
            bgRef = self.getBackgroundReference(patchRef)
            if not bgRef:
                return None
        tempExpName = self.getTempExpDatasetName()
        for tempExpRef in refList:
            if not tempExpRef.datasetExists(tempExpName):
                self.log.warn("Could not find %s %s; skipping it" % (tempExpName, tempExpRef.dataId))
                continue

            tempExp = tempExpRef.get(tempExpName, immediate=True)
            maskedImage = tempExp.getMaskedImage()
            imageScaler = self.scaleZeroPoint.computeImageScaler(exposure=tempExp, dataRef=tempExpRef)
            try:
                imageScaler.scaleMaskedImage(maskedImage)
            except Exception, e:
                self.log.warn("Scaling failed for %s (skipping it): %s" % (tempExpRef.dataId, e))
                continue
            statObj = afwMath.makeStatistics(maskedImage.getVariance(), maskedImage.getMask(),
                afwMath.MEANCLIP, statsCtrl)
            meanVar, meanVarErr = statObj.getResult(afwMath.MEANCLIP);
            weight = 1.0 / float(meanVar)
            if not numpy.isfinite(weight):
                self.log.warn("Non-finite weight for %s: skipping" % (tempExpRef.dataId,))
                continue
            self.log.info("Weight of %s %s = %0.3f" % (tempExpName, tempExpRef.dataId, weight))

            bgModel = None
            if self.config.doMatchBackgrounds:
                bgModel = self.matchBackgrounds.run(bgRef, tempExp, tract.getBBox(), inPlace=False)

                if self.debug:
                    suffix = "%s-%s" % (tempExpRef.dataId['visit'], patchRef.dataId['patch'])
                    tempExp.writeFits("orig-%s.fits" % suffix) # The warp, as it comes to us

                    bgImage = bgModel.getImageF(maskedImage.getBBox(afwImage.PARENT))
                    bgImage.writeFits("bgModel-%s.fits" % suffix) # The expanded background model
                    bgModel.getStatsImage().writeFits("bgStats-%s.fits" % suffix) # Background stats image

                    diffImage = bgRef.clone()
                    diffImage.getMaskedImage().__isub__(tempExp.getMaskedImage())
                    diffImage.writeFits("diff-%s.fits" % suffix) # bgRef - warp

                    warpImage = tempExp.clone()
                    warpImage.getMaskedImage().__iadd__(bgImage)
                    warpImage.writeFits("warp-%s.fits" % suffix) # This is what gets coadded

                    warpImage.getMaskedImage().__isub__(bgRef.getMaskedImage())
                    warpImage.writeFits("check-%s.fits" % suffix) # This should be zero

                    bgRef.writeFits("bgRef-%s.fits" % suffix)

                    del bgImage, diffImage, warpImage

            del maskedImage
            del tempExp

            tempExpRefList.append(tempExpRef)
            weightList.append(weight)
            imageScalerList.append(imageScaler)
            bgModelList.append(bgModel)

        return Struct(tempExpRefList=tempExpRefList, weightList=weightList,
                      imageScalerList=imageScalerList, bgModelList=bgModelList)

    def getBackgroundReference(self, patchRef):
        name = self.config.coaddName + "Coadd_bgRef"
        if not patchRef.datasetExists(name):
            return None
        bgRef = patchRef.get(name, immediate=True)
        model = Background.fromImage(self.config.background, bgRef.getMaskedImage())
        image = model.getImage()
        refImage = bgRef.getMaskedImage().getImage()
        refMask = bgRef.getMaskedImage().getMask()
        refMask.getArray()[:] = numpy.where(numpy.isnan(refImage.getArray()), 0, refMask.getArray())
        refImage.getArray()[:] = numpy.where(numpy.isnan(refImage.getArray()), image.getArray(), refImage.getArray())
        return bgRef

    def assemble(self, skyInfo, tempExpRefList, imageScalerList, weightList, bgModelList=None):
        """Assemble a coadd from input warps

        The assembly is performed over small areas on the image at a time, to
        conserve memory usage.

        @param skyInfo: Patch geometry information, from getSkyInfo
        @param tempExpRefList: List of data references to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching, or None
        @return coadded exposure
        """
        tempExpName = self.getTempExpDatasetName()
        self.log.info("Assembling %d %s for tract %s patch %s" %
                      (len(tempExpRefList), tempExpName, skyInfo.tractInfo.getId(),
                       skyInfo.patchInfo.getIndex()))

        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(self.getBadPixelMask())
        statsCtrl.setNanSafe(True)
        statsCtrl.setCalcErrorFromInputVariance(True)

        if self.config.doSigmaClip:
            statsFlags = afwMath.MEANCLIP
        else:
            statsFlags = afwMath.MEAN

        if bgModelList is None:
            bgModelList = [None]*len(tempExpRefList)

        coaddExposure = afwImage.ExposureF(skyInfo.bbox, skyInfo.wcs)
        coaddExposure.setCalib(self.scaleZeroPoint.getCalib())
        coaddExposure.getInfo().setCoaddInputs(self.inputRecorder.makeCoaddInputs())
        self.assembleMetadata(coaddExposure, tempExpRefList, weightList)
        coaddMaskedImage = coaddExposure.getMaskedImage()
        subregionSizeArr = self.config.subregionSize
        subregionSize = afwGeom.Extent2I(subregionSizeArr[0], subregionSizeArr[1])
        for subBBox in _subBBoxIter(skyInfo.bbox, subregionSize):
            try:
                self.assembleSubregion(coaddExposure, subBBox, tempExpRefList, imageScalerList,
                                       weightList, bgModelList, statsFlags, statsCtrl)
            except Exception, e:
                self.log.fatal("Cannot compute coadd %s: %s" % (subBBox, e,))
        self.log.info("Done coadding")
        coaddUtils.setCoaddEdgeBits(coaddMaskedImage.getMask(), coaddMaskedImage.getVariance())

        return coaddExposure

    def assembleMetadata(self, coaddExposure, tempExpRefList, weightList):
        """Set the metadata for the coadd

        This basic implementation simply sets the filter from the
        first input.

        @param coaddExposure: The target image for the coadd
        @param tempExpRefList: List of data references to tempExp
        @param weightList: List of weights
        """
        tempExpName = self.getTempExpDatasetName()
        # We load a single pixel of each coaddTempExp, because we just want to get at the metadata
        # (and we need more than just the PropertySet that contains the header).  See #2777.
        bbox = afwGeom.Box2I(afwGeom.Point2I(0,0), afwGeom.Extent2I(1,1))
        first = True
        coaddInputs = coaddExposure.getInfo().getCoaddInputs()
        for tempExpRef, weight in zip(tempExpRefList, weightList):
            tempExp = tempExpRef.get(tempExpName + "_sub", bbox=bbox, imageOrigin="LOCAL", immediate=True)
            if first:
                coaddExposure.setFilter(tempExp.getFilter())
                first = False
            self.inputRecorder.addVisitToCoadd(coaddInputs, tempExp, weight)
        coaddInputs.visits.sort()
        if self.config.doPsfMatch:
            psf = self.config.modelPsf.apply(coaddExposure.getWcs())
        else:
            psf = CoaddPsf(coaddInputs.ccds, coaddExposure.getWcs())
        coaddExposure.setPsf(psf)

    def assembleSubregion(self, coaddExposure, bbox, tempExpRefList, imageScalerList, weightList,
                          bgModelList, statsFlags, statsCtrl):
        """Assemble the coadd for a sub-region

        @param coaddExposure: The target image for the coadd
        @param bbox: Sub-region to coadd
        @param tempExpRefList: List of data reference to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching
        @param statsFlags: Statistic for coadd
        @param statsCtrl: Statistics control object for coadd
        """
        self.log.info("Computing coadd over %s" % bbox)
        tempExpName = self.getTempExpDatasetName()
        coaddMaskedImage = coaddExposure.getMaskedImage()
        coaddView = afwImage.MaskedImageF(coaddMaskedImage, bbox, afwImage.PARENT, False)
        maskedImageList = afwImage.vectorMaskedImageF() # [] is rejected by afwMath.statisticsStack
        for tempExpRef, imageScaler, bgModel in zip(tempExpRefList, imageScalerList, bgModelList):
            exposure = tempExpRef.get(tempExpName + "_sub", bbox=bbox, imageOrigin="PARENT")
            maskedImage = exposure.getMaskedImage()
            imageScaler.scaleMaskedImage(maskedImage)

            if self.config.doMatchBackgrounds:
                self.log.info("Applying background matching model to %s" % (tempExpRef.dataId,))
                self.matchBackgrounds.apply(maskedImage, bgModel, bbox, inPlace=True)

            maskedImageList.append(maskedImage)

        with self.timer("stack"):
            self.log.info("Stacking %d inputs" % len(maskedImageList))
            coaddSubregion = afwMath.statisticsStack(maskedImageList, statsFlags, statsCtrl, weightList)
            self.log.info("Done stacking")

        coaddView <<= coaddSubregion
        self.log.info("Copied in")

# XXX Old stuff, not used because of ParallelAssembleCoaddTask
#    def subtractBackground(self, exp):
#        # XXX these results from individual patches need to be merged
#        bgModel = getBackground(exp.getMaskedImage(), self.config.background)
#        bgImage = bgModel.getImageF(exp.getMaskedImage().getBBox(afwImage.PARENT),
#                                    self.config.background.algorithm, self.config.background.undersampleStyle)
#        exp.getMaskedImage().__isub__(bgImage)
#        return bgModel

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", cls.ConfigClass().coaddName + "Coadd_tempExp",
                               help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=AssembleCoaddDataIdContainer)
        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectDataIdContainer)
        return parser

class ParallelAssembleCoaddTask(SimpleAssembleCoaddTask):
    """Assemble a coadd from a set of coaddTempExp
    """

    def run(self, patchRefList, selectedDataList=[], visitKeys=["visit"], coaddName="deep"):


        # XXX debugging
        import itertools
        goodPatches = set(("%d,%d" % xy for xy in itertools.product(range(5,8),range(5,8))))
        patchRefList = [patchRef for patchRef in patchRefList if patchRef.dataId["patch"] in goodPatches]


        skymap = patchRefList[0].get(coaddName + "Coadd_skyMap")
        tractIdSet = set([patchRef.dataId["tract"] for patchRef in patchRefList])
        assert len(tractIdSet) == 1
        tract = skymap[tractIdSet.pop()]

        patchPool = Pool("assembleCoadd-patch")
        patchPool.storeSet(visitKeys=visitKeys, coaddName=coaddName)

        patchDataListList = patchPool.map(self.preparePatch, patchRefList, selectedDataList)

        # Pull out data for each visit
        visitData = {}
        for patchDataList, patchRef in zip(patchDataListList, patchRefList):
            if patchDataList is None:
                continue
            for patchData in patchDataList:
                visit = patchData.visit
                if visit in visitData:
                    visitData[visit].append(patchData)
                else:
                    visitData[visit] = [patchData]

        visitPool = Pool("assembleCoadd-visit") # Pool for processing visits
        visitPool.cacheClear()
        visitDataList = visitPool.map(self.mergeVisit, visitData.values(), tract)
        visitList = visitData.keys()
        if self.debug:
            for visit, data in zip(visitList, visitDataList):
                visitName = "+".join(map(str, visit))
                data.bgModel.getStatsImage().writeFits("bgMerge-%s.fits" % visitName)
                pickle.dump(data.bgModel, open("bgMerge-%s.pickle" % visitName, "w"))

        bgModelDict = dict(zip(visitList, [v.bgModel for v in visitDataList]))
        weightDict = dict(zip(visitList, [v.weight for v in visitDataList]))
        self.log.info("Visit weights: %s" % (weightDict,))
        del visitPool

        bgModelList = patchPool.mapToPrevious(self.assemblePatch, patchRefList, bgModelDict, weightDict)

        if self.config.doBackgroundSubtraction:
            bgModelList = [bgModel for bgModel in bgModelList if bgModel is not None]
            bgModel = reduce(lambda x,y: x.merge(y), bgModelList[1:], bgModelList[0])
            if self.debug:
                pickle.dump(bgModel, open("bg.pickle", "w"))
                bgModel.getStatsImage().writeFits("bg.fits")
        else:
            bgModel = None

        patchPool.mapToPrevious(self.completePatch, patchRefList, bgModel)

        return patchPool # Has the coadd in memory

    def preparePatch(self, cache, patchRef, selectDataList):
        cache.inputData = None
        skyInfo = self.getSkyInfo(patchRef)
        calExpRefList = self.selectExposures(patchRef, skyInfo, selectDataList=selectDataList)




        # XXXX dirty debugging hack
#        calExpRefList = [ref for ref in calExpRefList if ref.dataId['visit'] in (904014, 904010)]





        if len(calExpRefList) == 0:
            self.log.warn("No exposures to coadd for %s" % (patchRef.dataId,))
            return
        self.log.info("Coadding %d exposures for %s" % (len(calExpRefList), patchRef.dataId))

        tempExpRefList = self.getTempExpRefList(patchRef, calExpRefList)
        inputData = self.prepareInputs(patchRef, tempExpRefList, skyInfo.tractInfo)
        if not inputData:
            return None
        tempExpRefList = inputData.tempExpRefList
        self.log.info("Found %d %s for %s" % (len(inputData.tempExpRefList), self.getTempExpDatasetName(),
                                              patchRef.dataId))
        if len(inputData.tempExpRefList) == 0:
            self.log.warn("No coadd temporary exposures found for %s" % (patchRef.dataId,))
            return

        cache.calExpRefList = calExpRefList
        cache.inputData = inputData
        cache.skyInfo = skyInfo

        for ref, bgModel in zip(inputData.tempExpRefList, inputData.bgModelList):
            visit = tuple((ref.dataId[k] for k in cache.visitKeys))
            if self.debug:
                bgModel.getStatsImage().writeFits("bgMeas-%s-%s.fits" %
                                                  ("+".join(map(str, visit)),
                                                   ",".join(map(str, skyInfo.patchInfo.getIndex()))))

        return [Struct(visit=tuple((ref.dataId[k] for k in cache.visitKeys)), dataId=ref.dataId,
                       bgModel=bgModel, patch=skyInfo.patchInfo, weight=weight) for
                ref, bgModel, weight in zip(inputData.tempExpRefList, inputData.bgModelList,
                                            inputData.weightList)]

    def mergeVisit(self, cache, patchDataList, tract):
        patchList = [data.patch for data in patchDataList]
        bgModelList = [data.bgModel for data in patchDataList]
        weights = numpy.array([data.weight for data in patchDataList])
        bgModel = self.matchBackgrounds.merge(bgModelList) if self.config.doMatchBackgrounds else None
        return Struct(bgModel=bgModel,
                      weight=weights.mean(), # XXX weighted mean?
                      )

    def assemblePatch(self, cache, patchRef, bgModelDict, weightDict):
        cache.coaddExp = None
        if not cache.inputData:
            return None
        visitList = [tuple((ref.dataId[k] for k in cache.visitKeys)) for ref in cache.inputData.tempExpRefList]
        bgModelList = [bgModelDict[visit] for visit in visitList]
        weightList = [weightDict[visit] for visit in visitList]

        if self.debug:
            for tempExpRef, imageScaler, bgModel in zip(cache.inputData.tempExpRefList,
                                                        cache.inputData.imageScalerList,
                                                        bgModelList):
                tempExp = tempExpRef.get(cache.coaddName + "Coadd_tempExp", immediate=True)
                tempImage = tempExp.getMaskedImage()
                imageScaler.scaleMaskedImage(tempImage)
                if self.config.doMatchBackgrounds:
                    self.matchBackgrounds.apply(tempImage, bgModel, tempImage.getBBox(afwImage.PARENT),
                                                inPlace=True)
                visit = tuple((ref.dataId[k] for k in cache.visitKeys))
                patchIndex = cache.skyInfo.patchInfo.getIndex()
                tempExp.writeFits("input-%s-%s.fits" % ("+".join(map(str, visit)),
                                                        ",".join(map(str, patchIndex))))
                del tempExp, tempImage, visit, patchIndex

        coaddExp = self.assemble(cache.skyInfo, cache.inputData.tempExpRefList,
                                 cache.inputData.imageScalerList, weightList, bgModelList)

        if self.config.doInterp:
            self.log.info("Interpolating over bad pixels")
            fwhmPixels = self.config.interpFwhm / cache.skyInfo.wcs.pixelScale().asArcseconds()
            self.interpImage.interpolateOnePlane(maskedImage=coaddExp.getMaskedImage(), planeName="EDGE",
                                                 fwhmPixels=fwhmPixels)

        if self.config.doBackgroundSubtraction:
            # XXX Should be a two-step process, with a detection stage in-between;
            # simply running detection might work.
            self.log.info("Measuring coadd background")
            bgModel = Background.fromImage(self.config.background, coaddExp,
                                           box=cache.skyInfo.tractInfo.getBBox())
            if self.debug:
                pickle.dump(bgModel, open("bgSingle-%s.pickle" % (patchRef.dataId['patch'],), "w"))
                bgModel.getStatsImage().writeFits("bgSingle-%s.fits" % (patchRef.dataId['patch'],))
                bgImage = bgModel.getImageF(coaddExp.getMaskedImage().getBBox(afwImage.PARENT))
                bgImage.writeFits("bgImage-%s.fits" % (patchRef.dataId['patch'],))
        else:
            bgModel = None

        cache.coaddExp = coaddExp

        if self.debug:
            coaddExp.writeFits("coadd-%s.fits" % (patchRef.dataId['patch'],))

        return bgModel

    def completePatch(self, cache, patchRef, bgModel):
        if not cache.coaddExp:
            return
        if bgModel:
            bgImage = bgModel.getImageF(cache.skyInfo.bbox)
            if self.debug:
                bgImage.writeFits("bgApply-%s.fits" % (patchRef.dataId['patch'],))
            cache.coaddExp.getMaskedImage().__isub__(bgImage)
            if self.config.doWrite:
                self.writeCoaddOutput(patchRef, bgImage, "bg")

        self.writeCoaddOutput(patchRef, cache.coaddExp)


class StackConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    visitKeys = ListField(dtype=str, default=["visit"], doc="dataId keys that identify an exposure")
    select = ConfigurableField(target=BaseSelectImagesTask, doc="Select images to process")
    doWarp = Field(dtype=bool, default=True, doc="Warp images?")
    makeCoaddTempExp = ConfigurableField(target=MakeCoaddTempExpTask, doc="Warp images to sky")
    doBackgroundReference = Field(dtype=bool, default=True, doc="Build background reference?")
    backgroundReference = ConfigurableField(target=BackgroundReferenceTask, doc="Build background reference")
    assembleCoadd = ConfigurableField(target=ParallelAssembleCoaddTask, doc="Assemble warps into coadd")
    doProcessCoadd = Field(dtype=bool, default=True, doc="Process resulting coadd?")
    processCoadd = ConfigurableField(target=SubaruProcessCoaddTask, doc="Detection and measurement on coadd")
    doOverwriteCoadd = Field(dtype=bool, default=False, doc="Overwrite coadd?")
    doOverwriteOutput = Field(dtype=bool, default=False, doc="Overwrite processing outputs?")

    def setDefaults(self):
        self.select.retarget(WcsSelectImagesTask)
        self.makeCoaddTempExp.select.retarget(NullSelectImagesTask)
        self.backgroundReference.select.retarget(NullSelectImagesTask)
        self.assembleCoadd.select.retarget(NullSelectImagesTask)
        self.makeCoaddTempExp.doOverwrite = False
        self.processCoadd.detection.thresholdType = "pixel_stdev"
        self.assembleCoadd.doMatchBackgrounds = True
        self.makeCoaddTempExp.bgSubtracted = False

    def validate(self):
        if self.makeCoaddTempExp.coaddName != self.coaddName:
            raise RuntimeError("makeCoaddTempExp.coaddName and coaddName don't match")
        if self.assembleCoadd.coaddName != self.coaddName:
            raise RuntimeError("assembleCoadd.coaddName and coaddName don't match")
        if self.backgroundReference.coaddName != self.coaddName:
            raise RuntimeError("backgroundReference.coaddName and coaddName don't match")

class TractDataIdContainer(CoaddDataIdContainer):
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList

        It's difficult to make a data reference that merely points to an entire
        tract: there is no data product solely at the tract level.  Instead, we
        generate a list of data references for patches within the tract.
        """
        datasetType = namespace.config.coaddName + "Coadd"
        validKeys = set(["tract", "filter"])

        getPatchRefList = lambda tract: [namespace.butler.dataRef(datasetType=datasetType, tract=tract.getId(),
                                                                  filter=dataId["filter"],
                                                                  patch="%d,%d" % patch.getIndex()) for
                                         patch in tract]

        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract"):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(None, "--id must include " + key)

            skymap = self.getSkymap(namespace, datasetType)

            # tract is required; iterate over it if not provided
            if not "tract" in dataId:
                addList = [getPatchRefList(tract) for tract in skymap]
            else:
                addList = [getPatchRefList(skymap[dataId["tract"]])]
            self.refList += addList

class StackTaskRunner(CoaddTaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        if DEBUGGING:
            parsedCmd.selectId.dataList = None # XXX dirty debugging hack
        return CoaddTaskRunner.getTargetList(parsedCmd, butler=parsedCmd.butler, **kwargs)

class StackTask(PbsPoolTask):
    ConfigClass = StackConfig
    _DefaultName = "stacker" # "stack" conflicts with hscMosaic's StackTask.
    RunnerClass = StackTaskRunner

    def __init__(self, *args, **kwargs):
        super(StackTask, self).__init__(*args, **kwargs)
        self.makeSubtask("select")
        self.makeSubtask("makeCoaddTempExp")
        self.makeSubtask("backgroundReference")
        self.makeSubtask("assembleCoadd")
        self.makeSubtask("processCoadd")

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
        """
        Selection references are not cheap (reads Wcs), so are generated
        only if we're not doing a PBS submission.
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        # We don't want to be reading all the WCSes if we're only in the act of submitting to PBS
        SelectContainerClass = DataIdContainer if doPbs else SelectDataIdContainer
        if DEBUGGING:
            SelectContainerClass = DataIdContainer ### XXX dirty debugging hack
        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectContainerClass)
        return parser

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = len(parsedCmd.selectId.refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, patchRefList, butler, selectDataList=[]):
        """Run stacking on a tract

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param tractRef: Data reference for tract
        @param selectDataList: List of SelectStruct for inputs
        """
        import cPickle as pickle
        if DEBUGGING:
            self.log.warn("Reading pickled data...")
            patchRefList, selectDataList, coaddData = pickle.load(open("stack.pickle"))
        else:
            pool = Pool("stacker")
            pool.storeSet(warpType=self.config.coaddName + "Coadd_tempExp",
                          coaddType=self.config.coaddName + "Coadd")
            if self.config.doWarp:
                warpData = [Struct(patchRef=patchRef, selectDataList=selectDataList) for
                            patchRef in patchRefList]
                selectedData = pool.map(self.warp, warpData)
            else:
                # Faster not to do this in parallel, avoiding the overhead
                selectedData = [self.selectExposures(patchRef, selectDataList) for patchRef in patchRefList]
            if self.config.doBackgroundReference:
                self.backgroundReference.run(patchRefList, selectDataList)

            refNamer = lambda patchRef: tuple(map(int, patchRef.dataId["patch"].split(",")))
            lookup = dict(zip(map(refNamer, patchRefList), selectedData))
            coaddData = [Struct(patchRef=patchRef, selectDataList=lookup[refNamer(patchRef)]) for
                         patchRef in patchRefList]

            pickle.dump((patchRefList, selectDataList, coaddData), open("stack.pickle", "w"))

        patchPool = self.assembleCoadd.run(patchRefList, selectDataList, visitKeys=self.config.visitKeys,
                                           coaddName=self.config.coaddName)

        patchPool.storeSet(warpType=self.config.coaddName + "Coadd_tempExp",
                           coaddType=self.config.coaddName + "Coadd")
        if self.config.doProcessCoadd:
            patchPool.mapToPrevious(self.process, patchRefList)

    def warp(self, cache, data):
        """Warp all images for a patch

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param selectDataList: List of SelectStruct for inputs
        @return selectDataList with non-overlapping elements removed
        """
        patchRef = data.patchRef
        selectDataList = data.selectDataList
        self.log.info("%s: Start warping %s" % (NODE, patchRef.dataId))
        selectDataList = self.selectExposures(patchRef, selectDataList)
        self.makeCoaddTempExp.run(patchRef, selectDataList)
        self.log.info("%s: Finished warping %s" % (NODE, patchRef.dataId))
        return selectDataList

    def selectExposures(self, patchRef, selectDataList):
        """Select exposures to operate upon, via the SelectImagesTask

        This is very similar to CoaddBaseTask.selectExposures, except we return
        a list of SelectStruct (same as the input), so we can plug the results into
        future uses of SelectImagesTask.
        """
        key = lambda dataRef: tuple(dataRef.dataId[k] for k in sorted(dataRef.dataId.keys()))
        inputs = dict((key(select.dataRef), select) for select in selectDataList)
        skyMap = patchRef.get(self.config.coaddName + "Coadd_skyMap")
        tract = skyMap[patchRef.dataId["tract"]]
        patch = tract[(tuple(int(i) for i in patchRef.dataId["patch"].split(",")))]
        bbox = patch.getOuterBBox()
        wcs = tract.getWcs()
        cornerPosList = afwGeom.Box2D(bbox).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        dataRefList = self.select.runDataRef(patchRef, coordList, selectDataList=selectDataList).dataRefList
        return [inputs[key(dataRef)] for dataRef in dataRefList]

    def process(self, cache, patchRef):
        """Process (detection and measurement) image for a patch

        Only slave nodes execute this method.
        """
        coadd = cache.coaddExp
        if coadd is None:
            self.log.info("%s: No coadd to process for %s" % (NODE, patchRef.dataId))
            return
        if (not self.config.doOverwriteOutput and patchRef.datasetExists(cache.coaddType + "_src") and
            patchRef.datasetExists(cache.coaddType + "_calexp")):
            self.log.info("%s: Coadd outputs exist for %s" % (NODE, patchRef.dataId))
            return

        self.log.info("%s: Start processing %s" % (NODE, patchRef.dataId))
        try:
            self.processCoadd.process(patchRef, coadd)
        except LsstCppException as e:
            self.log.warn("LsstCppException %s" % NODE)
            if (isinstance(e.message, InvalidParameterException) and
                re.search("St. dev. must be > 0:", e.message.what())):
                # All the good pixels are outside the area of interest
                self.log.warn("No usable area for detection: %s" % patchRef.dataId)
                return
            raise
        self.log.info("%s: Finished processing %s" % (NODE, patchRef.dataId))

    def writeMetadata(self, dataRef):
        pass


"""
StackLauncher:
* Inputs: coadd name (for SkyMap retrieval), PBS stuff, exposure list, reference exposure
* Calculate overlaps: for each input, determine overlapping tract,patch: {visit/ccd --> tract/patch}
* Invert overlaps: {tract/patch --> visit/ccd}
* Determine background reference exposures (for now, use whatever's provided in refId; this is OK for single fields, but eventually we need a means of selecting reference exposures over an extended area; some ideas at https://dev.lsstcorp.org/trac/ticket/2741)
* Persist overlaps and background reference determinations
 save {tract/patch --> ref visit/ccd}
* For each tract/patch, launch StackTask

StackTask:
* Inputs: coadd name, tract/patch, specific exposure list, single ref exposure
* For each input: makeCoaddTempExp
* assembleCoadd with background matching to ref exposure
* (Could do a detection/photometry on the coadd right now, but background subtraction not the best)
* Generate summary statistics required for BackgroundSubtractionTask

BackgroundSubtractionTask:
* Inputs: coadd name, tract, patch list, ref exposure
* Determine background over the patches of the tract that share a common ref exposure id
* For each patch: write background model

StackFit = processCoadd:
* Subtract background, detect and measure
* Measure same apertures for different filters

ResolveDupes = done automatically in forcedPhotCoadd, but processCoadd dupe resolution requires LSST feature (#2893):
* Inputs: All sources measured from coadds
* Resolve duplicates (patch overlaps, tract overlaps)

ForcedPhotometry = forcedPhotCoadd:
* Inputs: coadd name, --id
* Get reference source list, measure with forced apertures

MultiFit = coming....:
* Inputs: coadd name, tract/patch, exposure list
* Get reference source list
* For each source in patch: determine relevant exposures, measure simultaneously across exposures


"""
