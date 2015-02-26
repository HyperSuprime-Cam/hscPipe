import re
import argparse

import numpy

from lsst.geom import convexHull

import lsst.pex.exceptions as pexExceptions
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.afw.table as afwTable
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
from lsst.afw.fits.fitsLib import FitsError
import lsst.coadd.utils as coaddUtils
from lsst.meas.algorithms import CoaddPsf, makeCoaddApCorrMap
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pex.exceptions import LsstCppException, InvalidParameterException
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, SelectDataIdContainer, CoaddTaskRunner
from lsst.pipe.tasks.selectImages import BaseSelectImagesTask, WcsSelectImagesTask, BaseExposureInfo
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.tasks.assembleCoadd import (AssembleCoaddTask, AssembleCoaddConfig, _subBBoxIter,
                                           AssembleCoaddDataIdContainer,)
from hsc.pipe.tasks.processCoadd import SubaruProcessCoaddTask
from lsst.pipe.tasks.coaddHelpers import groupPatchExposures, getGroupDataRef
from lsst.pipe.base import Struct, DataIdContainer, ArgumentParser
from hsc.pipe.base.parallel import BatchPoolTask
from hsc.pipe.base.pool import Pool, abortOnError, NODE
from hsc.pipe.base.butler import getDataRef
from hsc.pipe.tasks.background import BackgroundReferenceTask, MatchBackgroundsTask
from lsst.meas.algorithms import SourceDetectionTask

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
                                 inputData.weightList, inputData.bgModelList, doClip=self.config.doSigmaClip)

        if self.config.doInterp:
            self.interpImage.interpolateOnePlane(coaddExp.getMaskedImage(), "NO_DATA", coaddExp.getPsf())

        if self.config.doWrite:
            self.writeCoaddOutput(dataRef, coaddExp)

        return Struct(coaddExposure=coaddExp)

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

    def prepareInputs(self, patchRef, refList):
        """Prepare the input warps for coaddition

        This involves measuring weights and constructing image scalers
        for each of the inputs.

        @param refList: List of data references to tempExp
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
        tempExpName = self.getTempExpDatasetName()
        for tempExpRef in refList:
            if not tempExpRef.datasetExists(tempExpName):
                self.log.warn("Could not find %s %s; skipping it" % (tempExpName, tempExpRef.dataId))
                continue

            tempExp = tempExpRef.get(tempExpName, immediate=True)
            maskedImage = tempExp.getMaskedImage()
            imageScaler = self.scaleZeroPoint.computeImageScaler(
                exposure = tempExp,
                dataRef = tempExpRef,
            )
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
                bgModel = self.matchBackgrounds.run(bgRef, tempExp, inPlace=False)

                if self.debug:
                    suffix = "%s-%s" % (tempExpRef.dataId['visit'], patchRef.dataId['patch'])
                    tempExp.writeFits("orig-%s.fits" % suffix)

                    bgImage = bgModel.getImageF(self.matchBackgrounds.config.interpolation,
                                                self.matchBackgrounds.config.undersampleStyle)
                    bgImage.writeFits("bgModel-%s.fits" % suffix)

                    diffImage = bgRef.clone()
                    diffImage.getMaskedImage().__isub__(tempExp.getMaskedImage())
                    diffImage.writeFits("diff-%s.fits" % suffix)

                    warpImage = tempExp.clone()
                    warpImage.getMaskedImage().__iadd__(bgImage)
                    warpImage.writeFits("warp-%s.fits" % suffix)

                    warpImage.getMaskedImage().__isub__(bgRef.getMaskedImage())
                    warpImage.writeFits("check-%s.fits" % suffix)

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
        return patchRef.get(self.config.coaddName + "Coadd_bgRef", immediate=True)

    def assemble(self, skyInfo, tempExpRefList, imageScalerList, weightList, bgModelList=None,
                 altMaskList=None, doClip=False, mask=None):
        """Assemble a coadd from input warps

        The assembly is performed over small areas on the image at a time, to
        conserve memory usage.

        @param skyInfo: Patch geometry information, from getSkyInfo
        @param tempExpRefList: List of data references to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching, or None
        @param altMaskList: List of alternate masks to use than those on disk, or None
        @param doClip: Use clipping when codding?
        @param mask: Mask to ignore when coadding
        @return coadded exposure
        """
        tempExpName = self.getTempExpDatasetName()
        self.log.info("Assembling %s %s" % (len(tempExpRefList), tempExpName))
        if mask is None:
            mask = self.getBadPixelMask()

        statsCtrl = afwMath.StatisticsControl()
        statsCtrl.setNumSigmaClip(self.config.sigmaClip)
        statsCtrl.setNumIter(self.config.clipIter)
        statsCtrl.setAndMask(mask)
        statsCtrl.setNanSafe(True)
        statsCtrl.setWeighted(True)
        statsCtrl.setCalcErrorFromInputVariance(True)
        for plane, threshold in self.config.maskPropagationThresholds.items():
            bit = afwImage.MaskU.getMaskPlane(plane)
            statsCtrl.setMaskPropagationThreshold(bit, threshold)

        if doClip:
            statsFlags = afwMath.MEANCLIP
        else:
            statsFlags = afwMath.MEAN

        if bgModelList is None:
            bgModelList = [None]*len(tempExpRefList)

        if altMaskList is None:
            altMaskList = [None]*len(tempExpRefList)

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
                                       weightList, bgModelList, altMaskList, statsFlags, statsCtrl)
            except Exception, e:
                self.log.fatal("Cannot compute coadd %s: %s" % (subBBox, e,))

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
        apCorrMap = makeCoaddApCorrMap(coaddInputs.ccds, coaddExposure.getBBox(afwImage.PARENT),
                                       coaddExposure.getWcs())
        coaddExposure.getInfo().setApCorrMap(apCorrMap)

    def assembleSubregion(self, coaddExposure, bbox, tempExpRefList, imageScalerList, weightList,
                          bgModelList, altMaskList, statsFlags, statsCtrl):
        """Assemble the coadd for a sub-region

        @param coaddExposure: The target image for the coadd
        @param bbox: Sub-region to coadd
        @param tempExpRefList: List of data reference to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching
        @param altMaskList: List alternate masks to use than on disk
        @param statsFlags: Statistic for coadd
        @param statsCtrl: Statistics control object for coadd
        """
        self.log.logdebug("Computing coadd over %s" % bbox)
        tempExpName = self.getTempExpDatasetName()
        coaddMaskedImage = coaddExposure.getMaskedImage()
        coaddView = afwImage.MaskedImageF(coaddMaskedImage, bbox, afwImage.PARENT, False)
        maskedImageList = afwImage.vectorMaskedImageF() # [] is rejected by afwMath.statisticsStack
        for tempExpRef, imageScaler, bgModel, altMask in zip(tempExpRefList, imageScalerList, bgModelList,
                                                              altMaskList):
            exposure = tempExpRef.get(tempExpName + "_sub", bbox=bbox, imageOrigin="PARENT")
            maskedImage = exposure.getMaskedImage()

            if altMask:
                altMaskSub = altMask.Factory(altMask, bbox, afwImage.PARENT)
                maskedImage.getMask().swap(altMaskSub)

            imageScaler.scaleMaskedImage(maskedImage)

            if self.config.doMatchBackgrounds:
                bgImage = bgModel.getImageF()
                bgImage.setXY0(coaddMaskedImage.getXY0())
                maskedImage += bgImage.Factory(bgImage, bbox, afwImage.PARENT, False)

            maskedImageList.append(maskedImage)

        with self.timer("stack"):
            coaddSubregion = afwMath.statisticsStack(
                maskedImageList, statsFlags, statsCtrl, weightList)

        coaddView <<= coaddSubregion

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", cls.ConfigClass().coaddName + "Coadd_tempExp",
                               help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=AssembleCoaddDataIdContainer)
        parser.add_id_argument("--selectId", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectDataIdContainer)
        return parser


def countMaskFromFootprint(mask,footprint, bitmask, ignoreMask):
    """Function to count the number of pixels with a specific mask in a footprint
    """
    bbox = footprint.getBBox()
    bbox.clip(mask.getBBox(afwImage.PARENT))
    fp = afwImage.MaskU(bbox)
    subMask = mask.Factory(mask, bbox, afwImage.PARENT)
    afwDet.setMaskFromFootprint(fp, footprint, bitmask)
    return numpy.logical_and((subMask.getArray() & fp.getArray()) > 0,
                             (subMask.getArray() & ignoreMask) == 0).sum()

class SafeClipAssembleCoaddConfig(SimpleAssembleCoaddConfig):
    clipDetection = ConfigurableField(target=SourceDetectionTask,
                                      doc="Detect sources on difference between unclipped and clipped coadd")
    minClipFootOverlap = Field(
        doc = "Minimum fractional overlap of clipped footprint with visit DETECTED to be clipped",
        dtype = float,
        default = 0.6
    )
    minClipFootOverlapSingle = Field(
        doc = "Minimum fractional overlap of clipped footprint with visit DETECTED to be " \
            "clipped when only one visit overlaps",
        dtype = float,
        default = 0.5
    )
    minClipFootOverlapDouble = Field(
        doc = "Minimum fractional overlap of clipped footprints with visit DETECTED to be " \
            "clipped when two visits overlap",
        dtype = float,
        default = 0.45
    )
    maxClipFootOverlapDouble = Field(
        doc = "Maximum fractional overlap of clipped footprints with visit DETECTED when " \
            "considering two visits",
        dtype = float,
        default = 0.15
    )
    minBigOverlap = Field(
        doc = "Minimum number of pixels in footprint to use DETECTED mask from the single visits " \
            "when labeling clipped footprints",
        dtype = int,
        default = 100
    )

    def setDefaults(self):
        Config.setDefaults(self)
        self.clipDetection.reEstimateBackground = False
        self.clipDetection.returnOriginalFootprints = False
        self.clipDetection.thresholdPolarity = "both"
        self.clipDetection.thresholdValue = 2
        self.clipDetection.nSigmaToGrow = 4
        self.clipDetection.minPixels = 4
        self.clipDetection.isotropicGrow = True
        self.clipDetection.thresholdType = "pixel_stdev"
        self.sigmaClip = 1.5
        self.clipIter = 3

class SafeClipAssembleCoaddTask(SimpleAssembleCoaddTask):
    """Assemble a coadd, trying to identify and remove objects/artifacts that appear
    only in a small number of visits
    """
    ConfigClass = SafeClipAssembleCoaddConfig

    def __init__(self, *args, **kwargs):
        super(SimpleAssembleCoaddTask, self).__init__(*args, **kwargs)
        schema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask("clipDetection", schema=schema)

    def assemble(self, skyInfo, tempExpRefList, imageScalerList, weightList, bgModelList,
                 altMaskList=None, doClip=False, mask=None):
        """Assemble the coadd for a region including clipping of artifacts

        Identify clipped regions by detecting objects on the difference between unclipped and clipped coadd
        and then flag these regions on the individual vists so they are ignored in the coaddtion process

        @param skyInfo: Patch geometry information, from getSkyInfo
        @param tempExpRefList: List of data reference to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching
        @param altMaskList: Not used here
        @param doClip: Not used here
        @param mask: Not used here
        return coadd exposure
        """
        exp = self.buildDifferenceImage(skyInfo, tempExpRefList, imageScalerList, weightList, bgModelList)
        mask = exp.getMaskedImage().getMask()
        mask.addMaskPlane("CLIPPED")
        maskClipValue = mask.getPlaneBitMask("CLIPPED")

        # Get clipped footprints
        result = self.detectClip(exp, tempExpRefList)
        self.log.info('Found %d clipped objects' % len(result.clipFootprints))

        # Create mask of the current clipped footprints
        maskClip = mask.Factory(mask.getBBox(afwImage.PARENT))
        afwDet.setMaskFromFootprintList(maskClip, result.clipFootprints, maskClipValue)

        # Go to individual visits for big footprints
        maskClipBig = maskClip.Factory(mask.getBBox(afwImage.PARENT))
        self.detectClipBig(result.tempExpClipList, result.clipFootprints, result.clipIndices, maskClipBig)

        maskClip |= maskClipBig

        # Now combine again using mean, ignoring CLIPPED pixels
        badMaskPlanes = self.config.badMaskPlanes[:]
        badMaskPlanes.append("CLIPPED")
        badPixelMask = afwImage.MaskU.getPlaneBitMask(badMaskPlanes)
        coaddExp = SimpleAssembleCoaddTask.assemble(self, skyInfo, tempExpRefList, imageScalerList,
                                                    weightList, bgModelList, result.tempExpClipList,
                                                    doClip=False, mask=badPixelMask)

        # Set the coadd CLIPPED mask from the footprints since currently pixels that are masked
        # do not get propagated
        maskExp = coaddExp.getMaskedImage().getMask()
        maskExp |= maskClip

        return coaddExp


    def buildDifferenceImage(self, skyInfo, tempExpRefList, imageScalerList, weightList,
                             bgModelList):
        """Return an exposure that contains the difference between and unclipped and clipped coadd

        @param assemble: Function used to create coadd
        @param skyInfo: Patch geometry information, from getSkyInfo
        @param tempExpRefList: List of data reference to tempExp
        @param imageScalerList: List of image scalers
        @param weightList: List of weights
        @param bgModelList: List of background models from background matching
        @return Difference image of unclipped and clipped coadd wrapped in an Exposure
        """
        coaddMean = SimpleAssembleCoaddTask.assemble(self, skyInfo, tempExpRefList, imageScalerList,
                                                     weightList, bgModelList, doClip=False)

        # Build the clipped coadd
        coaddClip = SimpleAssembleCoaddTask.assemble(self, skyInfo, tempExpRefList, imageScalerList,
                                                     weightList, bgModelList, doClip=True)

        coaddDiff = coaddMean.getMaskedImage().Factory(coaddMean.getMaskedImage())
        coaddDiff -= coaddClip.getMaskedImage()
        exp = afwImage.ExposureF(coaddDiff)
        exp.setPsf(coaddMean.getPsf())
        return exp


    def detectClip(self, exp, tempExpRefList):
        """Detect clipped regions on an exposure and set the mask on the individual tempExp masks

        @param exp: Exposure to run detection on
        @param tempExpRefList: List of data reference to tempExp
        @return struct containg:
        - clippedFootprints: list of clipped footprints
        - clippedIndices: indices for each clippedFootprint in tempExpRefList
        - tempExpClipList: list of new masks for tempExp
        """
        mask = exp.getMaskedImage().getMask()
        maskClipValue = mask.getPlaneBitMask("CLIPPED")
        maskDetValue = mask.getPlaneBitMask("DETECTED") | mask.getPlaneBitMask("DETECTED_NEGATIVE")
        fpSet = self.clipDetection.detectFootprints(exp, doSmooth=True, clearMask=True)
        # Merge positive and negative together footprints together
        fpSet.positive.merge(fpSet.negative)
        footprints = fpSet.positive
        self.log.info('Found %d potential clipped objects' % len(footprints.getFootprints()))
        ignoreMask = self.getBadPixelMask()

        clipFootprints = []
        clipIndices = []

        # build a list with a mask for each visit which can be modified with clipping information
        tempExpClipList = []
        for tmpExpRef in tempExpRefList:
            tmpExp = tmpExpRef.get(self.getTempExpDatasetName(), immediate=True)
            tmpExpMask = tmpExp.getMaskedImage().getMask()
            tempExpClipList.append(tmpExpMask)

        for footprint in footprints.getFootprints():
            nPixel = footprint.getArea()
            overlap = [] # hold the overlap with each visit
            maskList = [] # which visit mask match
            indexList = []# index of visit in global list
            for i,tmpExpMask in enumerate(tempExpClipList):
                # what is the overlap with the footprint
                ignore = countMaskFromFootprint(tmpExpMask, footprint, ignoreMask, 0x0)
                overlapDet = countMaskFromFootprint(tmpExpMask, footprint, maskDetValue,
                                                    ignoreMask)
                totPixel = nPixel - ignore

                # If we have more bad pixels than detection skip
                if ignore > overlapDet or totPixel == 0 or overlapDet == 0:
                    continue
                overlap.append(overlapDet/float(totPixel))
                maskList.append(tmpExpMask)
                indexList.append(i)

            overlap = numpy.array(overlap)
            if len(overlap) == 0:
                continue

            keep = False   # Should this footprint be marked as clipped?
            keepIndex = [] # Which tempExps does the clipped footprint belong to

            # If footprint only has one overlap use a lower threshold
            if len(overlap) == 1:
                if overlap[0] > self.config.minClipFootOverlapSingle:
                    keep = True
                    keepIndex = [0]
            else:
                # This is the general case for most clipped objects where only visit is responsible
                clipIndex = numpy.where(overlap > self.config.minClipFootOverlap)[0]
                if len(clipIndex) == 1:
                    keep=True
                    keepIndex = [clipIndex[0]]

                # Test if there are clipped objects that overlap two different visits
                clipIndex = numpy.where(overlap > self.config.minClipFootOverlapDouble)[0]
                if len(clipIndex) == 2 and len(overlap) >= 3:
                    clipIndexComp = numpy.where(overlap < self.config.minClipFootOverlapDouble)[0]
                    if numpy.max(overlap[clipIndexComp]) < self.config.maxClipFootOverlapDouble:
                        keep=True
                        keepIndex = clipIndex

            if not keep:
                continue

            for index in keepIndex:
                afwDet.setMaskFromFootprint(maskList[index], footprint, maskClipValue)

            clipIndices.append(numpy.array(indexList)[keepIndex])
            clipFootprints.append(footprint)

        return Struct(clipFootprints=clipFootprints, clipIndices=clipIndices,
                      tempExpClipList=tempExpClipList)


    def detectClipBig(self, tempExpClipList, clipFootprints, clipIndices, maskClipBig):
        """Find footprints from individual tempExp footprints for large footprints

        This will update maskClipBig and the tempExp masks

        @param tempExpClipList: List of tempExp mask's with clipping information
        @param clipFootprints: List of clipped footprints
        @param clipIndices: List of which entries in tempExpClipList each footprint belongs to
        @param maskClipValue: Mask of clipped objects from coadd
        @param maskClipValue: Mask value of clipped pixels
        @param maskClipValue: Mask value of detected pixels
        @return list of big footprints
        """
        ignoreMask = self.getBadPixelMask()
        maskClipValue = maskClipBig.getPlaneBitMask("CLIPPED")
        maskDetValue = maskClipBig.getPlaneBitMask("DETECTED")
        maskDetValue |= maskClipBig.getPlaneBitMask("DETECTED_NEGATIVE")
        for index,tmpExpMask in enumerate(tempExpClipList):

            # Create list of footprints from the DETECTED pixels
            maskVisitDet = tmpExpMask.Factory(tmpExpMask, tmpExpMask.getBBox(afwImage.PARENT),
                                              afwImage.PARENT, True)
            maskVisitDet &= maskDetValue
            visitFootprints = afwDet.FootprintSet(maskVisitDet, afwDet.Threshold(1))

            # build a mask of clipped footprints that are in this visit
            clippedFootprintsVisit = []
            for foot,clipIndex in zip(clipFootprints, clipIndices):
                if index not in clipIndex:
                    continue
                clippedFootprintsVisit.append(foot)
            maskVisitClip = maskVisitDet.Factory(maskVisitDet.getBBox(afwImage.PARENT))
            afwDet.setMaskFromFootprintList(maskVisitClip, clippedFootprintsVisit, maskClipValue)

            bigFootprintsVisit = []
            for foot in visitFootprints.getFootprints():
                if foot.getArea() < self.config.minBigOverlap:
                    continue
                nCount = countMaskFromFootprint(maskVisitClip, foot, maskClipValue, ignoreMask)
                # Does it overlap an existing big footprint?  For bright stars that can be split up
                # into multiple pieces we don't want to mask it out of every visit so we only use
                # the first one
                overlapExisting = countMaskFromFootprint(maskClipBig, foot, maskClipValue, 0x0)
                if nCount > self.config.minBigOverlap and overlapExisting == 0:
                    bigFootprintsVisit.append(foot)
                    afwDet.setMaskFromFootprint(maskClipBig, foot, maskClipValue)

            # Update single visit masks
            maskVisitClip.clearAllMaskPlanes()
            afwDet.setMaskFromFootprintList(maskVisitClip, bigFootprintsVisit, maskClipValue)
            tmpExpMask |= maskVisitClip



class StackConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    select = ConfigurableField(target=BaseSelectImagesTask, doc="Select images to process")
    makeCoaddTempExp = ConfigurableField(target=MakeCoaddTempExpTask, doc="Warp images to sky")
    doBackgroundReference = Field(dtype=bool, default=True, doc="Build background reference?")
    backgroundReference = ConfigurableField(target=BackgroundReferenceTask, doc="Build background reference")
    assembleCoadd = ConfigurableField(target=SafeClipAssembleCoaddTask, doc="Assemble warps into coadd")
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
        validKeys = set(["tract", "filter", "patch",])

        getPatchRefList = lambda tract: [namespace.butler.dataRef(datasetType=datasetType, tract=tract.getId(),
                                                                  filter=dataId["filter"],
                                                                  patch="%d,%d" % patch.getIndex()) for
                                         patch in tract]

        tractRefs = {} # Data references for each tract
        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract", "patch",):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(None, "--id must include " + key)

            skymap = self.getSkymap(namespace, datasetType)

            if "tract" in dataId:
                tractId = dataId["tract"]
                if tractId not in tractRefs:
                    tractRefs[tractId] = []
                if "patch" in dataId:
                    tractRefs[tractId].append(namespace.butler.dataRef(datasetType=datasetType, tract=tractId,
                                                                       filter=dataId['filter'],
                                                                       patch=dataId['patch']))
                else:
                    tractRefs[tractId] += getPatchRefList(skymap[tractId])
            else:
                tractRefs = dict((tract.getId(), tractRefs.get(tract.getId(), []) + getPatchRefList(tract)) for
                                 tract in skymap)

        self.refList = tractRefs.values()

class StackTaskRunner(CoaddTaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        kwargs["butler"] = parsedCmd.butler
        kwargs["selectIdList"] = [ref.dataId for ref in parsedCmd.selectId.refList]
        return [(parsedCmd.id.refList, kwargs),]

class StackTask(BatchPoolTask):
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
    def _makeArgumentParser(cls, doBatch=False, **kwargs):
        """
        Selection references are not cheap (reads Wcs), so are generated
        only if we're not doing a batch submission.
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        parser.add_id_argument("--selectId", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9")
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = len(parsedCmd.selectId.refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, tractPatchRefList, butler, selectIdList=[]):
        """Determine which tracts are non-empty before processing"""
        pool = Pool("tracts")
        pool.storeSet(butler=butler, skymap=butler.get(self.config.coaddName + "Coadd_skyMap"))
        tractIdList = []
        for patchRefList in tractPatchRefList:
            tractSet = set([patchRef.dataId["tract"] for patchRef in patchRefList])
            assert len(tractSet) == 1
            tractIdList.append(tractSet.pop())

        selectDataList = [data for data in pool.mapNoBalance(self.readSelection, selectIdList) if
                          data is not None]

        nonEmptyList = pool.mapNoBalance(self.checkTract, tractIdList, selectDataList)
        tractPatchRefList = [patchRefList for patchRefList, nonEmpty in
                             zip(tractPatchRefList, nonEmptyList) if nonEmpty]
        self.log.info("Non-empty tracts (%d): %s" % (len(tractPatchRefList),
                                                     [patchRefList[0].dataId["tract"] for patchRefList in
                                                      tractPatchRefList]))

        # Install the dataRef in the selectDataList
        for data in selectDataList:
            data.dataRef = getDataRef(butler, data.dataId, "calexp")

        # Process the non-empty tracts
        return [self.runTract(patchRefList, butler, selectDataList) for patchRefList in tractPatchRefList]

    @abortOnError
    def runTract(self, patchRefList, butler, selectDataList=[]):
        """Run stacking on a tract

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param tractRef: Data reference for tract
        @param selectDataList: List of SelectStruct for inputs
        """
        pool = Pool("stacker")
        pool.cacheClear()
        pool.storeSet(butler=butler, warpType=self.config.coaddName + "Coadd_tempExp",
                      coaddType=self.config.coaddName + "Coadd")
        patchIdList = [patchRef.dataId for patchRef in patchRefList]

        selectedData = pool.map(self.warp, patchIdList, selectDataList)
        if self.config.doBackgroundReference:
            self.backgroundReference.run(patchRefList, selectDataList)

        refNamer = lambda patchRef: tuple(map(int, patchRef.dataId["patch"].split(",")))
        lookup = dict(zip(map(refNamer, patchRefList), selectedData))
        coaddData = [Struct(patchId=patchRef.dataId, selectDataList=lookup[refNamer(patchRef)]) for
                     patchRef in patchRefList]
        pool.map(self.coadd, coaddData)

    def readSelection(self, cache, selectId):
        """Read Wcs of selected inputs

        This method only runs on slave nodes.

        This method is similar to SelectDataIdContainer.makeDataRefList,
        creating a Struct like a SelectStruct, except with a dataId instead
        of a dataRef (to ease MPI).

        @param cache: Pool cache
        @param selectId: Data identifier for selected input
        @return a SelectStruct with a dataId instead of dataRef
        """
        try:
            ref = getDataRef(cache.butler, selectId, "calexp")
            self.log.info("Reading Wcs from %s" % (selectId,))
            md = ref.get("calexp_md", immediate=True)
            wcs = afwImage.makeWcs(md)
            data = Struct(dataId=selectId, wcs=wcs, dims=(md.get("NAXIS1"), md.get("NAXIS2")))
        except pexExceptions.LsstCppException, e:
            if not isinstance(e.message, FitsError): # Unable to open file
                raise
            self.log.warn("Unable to construct Wcs from %s" % (selectId,))
            return None
        return data

    def checkTract(self, cache, tractId, selectIdList):
        """Check whether a tract has any overlapping inputs

        This method only runs on slave nodes.

        @param cache: Pool cache
        @param tractId: Data identifier for tract
        @param selectDataList: List of selection data
        @return whether tract has any overlapping inputs
        """
        skymap = cache.skymap
        tract = skymap[tractId]
        tractWcs = tract.getWcs()
        tractPoly = convexHull([tractWcs.pixelToSky(afwGeom.Point2D(coord)).getVector() for
                                coord in tract.getBBox().getCorners()])

        for selectData in selectIdList:
            if not hasattr(selectData, "poly"):
                wcs = selectData.wcs
                dims = selectData.dims
                box = afwGeom.Box2D(afwGeom.Point2D(0, 0), afwGeom.Point2D(*dims))
                selectData.poly = convexHull([wcs.pixelToSky(coord).getVector() for coord in box.getCorners()])
            if tractPoly.intersects(selectData.poly):
                return True
        return False

    def warp(self, cache, patchId, selectDataList):
        """Warp all images for a patch

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param selectDataList: List of SelectStruct for inputs
        @return selectDataList with non-overlapping elements removed
        """
        patchRef = getDataRef(cache.butler, patchId, cache.coaddType)
        selectDataList = self.selectExposures(patchRef, selectDataList)
        with self.logOperation("warping %s" % (patchRef.dataId,), catch=True):
            self.makeCoaddTempExp.run(patchRef, selectDataList)
        return selectDataList

    def coadd(self, cache, data):
        """Construct coadd for a patch and measure

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param selectDataList: List of SelectStruct for inputs
        """
        patchRef = getDataRef(cache.butler, data.patchId, cache.coaddType)
        selectDataList = data.selectDataList
        coadd = None
        with self.logOperation("coadding %s" % (patchRef.dataId,), catch=True):
            if self.config.doOverwriteCoadd or not patchRef.datasetExists(cache.coaddType):
                coaddResults = self.assembleCoadd.run(patchRef, selectDataList)
                if coaddResults is not None:
                    coadd = coaddResults.coaddExposure
            elif patchRef.datasetExists(cache.coaddType):
                self.log.info("%s: Reading coadd %s" % (NODE, patchRef.dataId))
                coadd = patchRef.get(cache.coaddType, immediate=True)

        if coadd is not None and (self.config.doOverwriteOutput or
                                  not patchRef.datasetExists(cache.coaddType + "_src") or
                                  not patchRef.datasetExists(cache.coaddType + "_calexp")):
                self.process(patchRef, coadd)

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

    def process(self, patchRef, coadd):
        with self.logOperation("processing %s" % (patchRef.dataId,), catch=True):
            try:
                self.processCoadd.run(patchRef)
            except LsstCppException as e:
                self.log.warn("LsstCppException %s" % NODE)
                if (isinstance(e.message, InvalidParameterException) and
                    re.search("St. dev. must be > 0:", e.message.what())):
                    # All the good pixels are outside the area of interest; allow to proceed
                    self.log.warn("No usable area for detection: %s" % patchRef.dataId)
                else:
                    raise

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
