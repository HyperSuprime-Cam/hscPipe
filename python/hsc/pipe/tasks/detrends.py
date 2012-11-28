#!/usr/bin/env python

import os
import numpy
import argparse

import pbasf2 as pbasf
from lsst.pex.config import Config, ConfigField, ConfigurableField, Field
from lsst.pipe.base import Task, Struct
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
import lsst.meas.algorithms as measAlg

import lsst.obs.subaru.isr as hscIsr

import hsc.pipe.base.butler as hscButler
from hsc.pipe.base.mpi import MpiTask, MpiArgumentParser, abortOnError, thisNode
from hsc.pipe.base.pbs import PbsArgumentParser, shCommandFromArgs

class DetrendStatsConfig(Config):
    """Parameters controlling background statistics"""
    stat = Field(doc="Statistic to use to estimate background (from lsst.afw.math)", dtype=int,
                   default=afwMath.MEANCLIP)
    clip = Field(doc="Clipping threshold for background", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for background", dtype=int, default=3)

class DetrendStatsTask(Task):
    """Measure statistic of the background"""
    ConfigClass = DetrendStatsConfig

    def run(self, exposureOrImage):
        """Measure a particular statistic on an image (of some sort).

        @param exposureOrImage    Exposure, MaskedImage or Image.
        @return Value of desired statistic
        """
        stats = afwMath.StatisticsControl(self.config.clip, self.config.iter,
                                          afwImage.MaskU.getPlaneBitMask("DETECTED"))
        try:
            image = exposureOrImage.getMaskedImage()
        except:
            try:
                image = exposureOrImage.getImage()
            except:
                image = exposureOrImage

        return afwMath.makeStatistics(image, self.config.stat, stats).getValue()


class DetrendCombineConfig(Config):
    rows = Field(doc="Number of rows to read at a time", dtype=int, default=128)
    maskDetected = Field(doc="Mask pixels about the detection threshold?", dtype=bool, default=True)
    combine = Field(doc="Statistic to use for combination (from lsst.afw.math)", dtype=int,
                    default=afwMath.MEANCLIP)
    clip = Field(doc="Clipping threshold for combination", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for combination", dtype=int, default=3)
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")

class DetrendCombineTask(Task):
    ConfigClass = DetrendCombineConfig

    def __init__(self, *args, **kwargs):
        super(DetrendCombineTask, self).__init__(*args, **kwargs)
        self.makeSubtask("stats")

    def run(self, sensorRefList, expScales=None, finalScale=None, outputName="postISRCCD"):
        width, height = self.getDimensions(sensorRefList)
        maskVal = afwImage.MaskU.getPlaneBitMask("DETECTED") if self.config.maskDetected else 0
        stats = afwMath.StatisticsControl(self.config.clip, self.config.iter, maskVal)

        # Combine images
        combined = afwImage.ImageF(width, height)
        numImages = len(sensorRefList)
        imageList = afwImage.vectorMaskedImageF(numImages)
        for start in range(0, height, self.config.rows):
            rows = min(self.config.rows, height - start)
            box = afwGeom.Box2I(afwGeom.Point2I(0, start), afwGeom.Extent2I(width, rows))
            subCombined = combined.Factory(combined, box)

            for i, sensorRef in enumerate(sensorRefList):
                exposure = sensorRef.get(outputName + "_sub", bbox=box)
                if expScales is not None:
                    self.applyScale(exposure, expScales[i])
                imageList[i] = exposure.getMaskedImage()

            self.combine(subCombined, imageList, stats)

        if finalScale is not None:
            background = self.stats.run(combined)
            self.log.log(self.log.INFO, "Measured background of stack is %f; adjusting to %f" %
                         (background, finalScale))
            combined *= finalScale / background

        return combined

    def getDimensions(self, sensorRefList, outputName="postISRCCD"):
        dimList = []
        for sensorRef in sensorRefList:
            md = sensorRef.get(outputName + "_md")
            dimList.append(afwGeom.Extent2I(md.get("NAXIS1"), md.get("NAXIS2")))
        return getSize(dimList)

    def applyScale(self, exposure, scale=None):
        if scale is not None:
            mi = exposure.getMaskedImage()
            mi /= scale

    def combine(self, target, imageList, stats):
        if False:
            # In-place stacks are now supported on LSST's afw, but not yet on HSC
            afwMath.statisticsStack(target, imageList, self.config.combine, stats)
        else:
            stack = afwMath.statisticsStack(imageList, self.config.combine, stats)
            target <<= stack.getImage()



def getCcdName(ccdId, ccdKeys):
    return tuple(ccdId[k] for k in ccdKeys)

def getCcdIdListFromExposures(expRefList, level="sensor"):
    """Determine a list of CCDs from exposure references

    This essentially inverts the exposure-level references (which
    provides a list of CCDs for each exposure), by providing
    a set of keywords that identify a CCD in the dataId, and a
    dataId list for each CCD.  Consider an input list of exposures
    [e1, e2, e3], and each exposure has CCDs c1 and c2.  Then this
    function returns:

        set(['ccd']),
        {(c1,): [e1c1, e2c1, e3c1], (c2,): [e1c2, e2c2, e3c2]}

    The latter part is a dict whose keys are tuples of the identifying
    values of a CCD (usually just the CCD number) and the values are
    lists of dataIds for that CCD in each exposure.  A missing dataId
    is given the value None.
    """
    expIdList = [[ccdRef.dataId for ccdRef in expRef.subItems(level)] for expRef in expRefList]

    # Determine what additional keys make a CCD from an exposure
    ccdKeys = set() # Set of keywords in the dataId that identify a CCD
    ccdNames = set() # Set of tuples which are values for each of the CCDs in an exposure
    for ccdIdList, expRef in zip(expIdList, expRefList):
        expKeys = set(expRef.dataId.keys())
        for ccdId in ccdIdList:
            keys = set(ccdId.keys()).difference(expKeys)
            if len(ccdKeys) == 0:
                ccdKeys = keys
            elif keys != ccdKeys:
                raise RuntimeError("Keys for CCD differ: %s vs %s" % (keys, ccdList.keys()))
            name = getCcdName(ccdId, ccdKeys)
            ccdNames.add(name)

    # Turn the list of CCDs for each exposure into a list of exposures for each CCD
    ccdLists = dict((k,[]) for k in ccdNames)
    for n, ccdIdList in enumerate(expIdList):
        for ccdId in ccdIdList:
            name = getCcdName(ccdId, ccdKeys)
            ccdLists[name].append(ccdId)
        for idList in ccdLists.values():
            if len(idList) == n:
                idList.append(None)

    return ccdKeys, ccdLists

class DetrendIdAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string):
        output = getattr(namespace, self.dest, {})
        for nameValue in values:
            name, sep, valueStr = nameValue.partition("=")
            if not valueStr:
                parser.error("%s value %s must be in form name=value" % (option_string, nameValue))
            output[name] = valueStr
        setattr(namespace, self.dest, output)

class DetrendArgumentParser(MpiArgumentParser):
    def __init__(self, calibName, *args, **kwargs):
        super(DetrendArgumentParser, self).__init__(*args, **kwargs)
        self.calibName = calibName
        self.add_argument("--detrendId", nargs="*", action=DetrendIdAction, default={},
                          help="identifiers for detrend, e.g., --detrendId version=1",
                          metavar="KEY=VALUE1[^VALUE2[^VALUE3...]")
    def parse_args(self, *args, **kwargs):
        namespace = super(DetrendArgumentParser, self).parse_args(*args, **kwargs)

        keys = namespace.butler.getKeys(self.calibName)
        parsed = {}
        for name, value in namespace.detrendId.items():
            if not name in keys:
                parser.error("%s is not a relevant detrend identifier key (%s)" % (name, keys))
            parsed[name] = keys[name](value)
        namespace.detrendId = parsed

        return namespace

class DetrendConfig(Config):
    isr = ConfigurableField(target=hscIsr.SubaruIsrTask, doc="ISR configuration")
    dateObs = Field(dtype=str, default="dateObs", doc="Key for observation date in exposure registry")
    dateCalib = Field(dtype=str, default="calibDate", doc="Key for detrend date in calib registry")
    filter = Field(dtype=str, default="filter", doc="Key for filter name in exposure/calib registries")
    combination = ConfigurableField(target=DetrendCombineTask, doc="Detrend combination configuration")
    def setDefaults(self):
        self.isr.doWrite = False

class DetrendTask(MpiTask):
    """Quite abstract (though not completely) base class for combining detrends"""
    ConfigClass = DetrendConfig

    def __init__(self, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(DetrendTask, self).__init__(**kwargs)
        self.makeSubtask("isr")
        self.makeSubtask("combination")

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        return DetrendArgumentParser(calibName=cls.calibName, name=cls._DefaultName,
                                     dataRefLevel="visit", *args, **kwargs)

    @abortOnError
    def runDataRefList(self, expRefList, doRaise=False):
        """All nodes execute this"""
        self.butler = self.parsedCmd.butler

        if self.rank == self.root:
            outputId = self.getOutputId(expRefList)
            ccdKeys, ccdIdLists = getCcdIdListFromExposures(expRefList, level="sensor")

            # Ensure we can generate filenames for each output
            for ccdName in ccdIdLists:
                dataId = dict(outputId.items() + [(k, ccdName[i]) for i, k in enumerate(ccdKeys)])
                try:
                    filename = self.butler.get(self.calibName + "_filename", dataId)
                except Exception, e:
                    raise RuntimeError("Unable to determine output filename from %s: %s" % (dataId, e))
        else:
            outputId = None
            ccdKeys, ccdIdLists = None, {}
        ccdIdLists = pbasf.Broadcast(self.comm, ccdIdLists, root=self.root)

        # Scatter: process CCDs independently
        data = self.scatterProcess(ccdKeys, ccdIdLists)

        # Gather: determine scalings
        if self.rank == self.root:
            scales = self.scale(ccdKeys, ccdIdLists, data)
        else:
            scales = None, None

        # Scatter: combine
        self.scatterCombine(outputId, ccdKeys, ccdIdLists, scales)

    def getOutputId(self, expRefList):
        expIdList = [expRef.dataId for expRef in expRefList]
        midTime = 0
        filterName = None
        for expId in expIdList:
            midTime += self.getMjd(expId)
            thisFilter = self.getFilter(expId)
            if filterName is None:
                filterName = thisFilter
            elif filterName != thisFilter:
                raise RuntimeError("Filter mismatch for %s: %s vs %s" % (expId, thisFilter, filterName))

        midTime /= len(expRefList)
        date = str(dafBase.DateTime(midTime, dafBase.DateTime.MJD).toPython().date())

        outputId = {self.config.filter: filterName, self.config.dateCalib: date}
        outputId.update(self.parsedCmd.detrendId)
        return outputId


    def getMjd(self, dataId):
        dateObs = dataId[self.config.dateObs]
        try:
            dt = dafBase.DateTime(dateObs)
        except:
            dt = dafBase.DateTime(dateObs + "T12:00:00.0Z")
        return dt.get(dafBase.DateTime.MJD)

    def getFilter(self, dataId):
        return dataId[self.config.filter]

    def scatterProcess(self, ccdKeys, ccdIdLists):
        """Scatter the processing among the nodes

        We scatter the data wider than the just the number of CCDs, to make
        full use of all available processors.  This necessitates piecing
        everything back together in the same format as ccdIdLists afterwards.
        """
        if self.rank == self.root:
            dataIdList = sum(ccdIdLists.values(), [])
        else:
            dataIdList = None

        self.log.info("Scatter processing on %s" % thisNode())
        resultList = pbasf.ScatterJob(self.comm, self.process, dataIdList, root=self.root)
        if self.rank == self.root:
            data = dict((ccdName, [None] * len(expList)) for ccdName, expList in ccdIdLists.items())
            indices = dict(sum([[(tuple(dataId.values()), (ccdName, expNum))
                                 for expNum, dataId in enumerate(expList)]
                                for ccdName, expList in ccdIdLists.items()], []))
            for dataId, result in zip(dataIdList, resultList):
                ccdName, expNum = indices[tuple(dataId.values())]
                data[ccdName][expNum] = result
        else:
            data = None
        return data

    def process(self, ccdId):
        self.log.info("Processing %s on %s" % (ccdId, thisNode()))
        sensorRef = hscButler.getDataRef(self.butler, ccdId)
        exposure = self.processSingle(sensorRef)
        self.processWrite(sensorRef, exposure)
        return self.processResult(exposure)

    def processSingle(self, dataRef):
        return self.isr.run(dataRef).exposure

    def processWrite(self, dataRef, exposure, outputName="postISRCCD"):
        dataRef.put(exposure, outputName)

    def processResult(self, exposure):
        return None

    def scale(self, ccdKeys, ccdIdLists, data):
        self.log.info("Scale on %s" % thisNode())
        return dict((name, Struct(ccdScale=None, expScales=[None] * len(ccdIdLists[name])))
                    for name in ccdIdLists.keys())

    def scatterCombine(self, outputId, ccdKeys, ccdIdLists, scales):
        self.log.info("Scatter combination on %s" % thisNode())
        if self.rank == self.root:
            data = [Struct(ccdIdList=ccdIdLists[ccdName], scales=scales[ccdName],
                           outputId=dict(outputId.items() + [(k,ccdName[i]) for i, k in enumerate(ccdKeys)]))
                    for ccdName in ccdIdLists.keys()]
        else:
            data = None
        pbasf.ScatterJob(self.comm, self.combine, data, root=self.root)

    def combine(self, struct):
        dataRefList = [hscButler.getDataRef(self.butler, dataId) for dataId in struct.ccdIdList]
        self.log.info("Combining %s on %s" % (struct.outputId, thisNode()))
        detrend = self.combination.run(dataRefList, expScales=struct.scales.expScales,
                                       finalScale=struct.scales.ccdScale)
        self.write(detrend, struct.outputId)

    def write(self, exposure, dataId):
        self.log.info("Writing %s on %s" % (dataId, thisNode()))
        self.butler.put(exposure, self.calibName, dataId)


class BiasConfig(DetrendConfig):
    pass

def biasOverrides(config):
    config.isr.doBias = False
    config.isr.doDark = False
    config.isr.doFlat = False
#   config.isr.doFringe = False


class BiasTask(DetrendTask):
    ConfigClass = BiasConfig
    _DefaultName = "bias"
    overrides = [biasOverrides]
    calibName = "bias"

class DarkConfig(DetrendConfig):
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe")

def darkOverrides(config):
    config.isr.doDark = False
    config.isr.doFlat = False
#   config.isr.doFringe = False

class DarkTask(BiasTask):
    ConfigClass = DarkConfig
    _DefaultName = "dark"
    overrides = [darkOverrides]
    calibName = "dark"

    def processSingle(self, sensorRef):
        exposure = super(DarkTask, self).processSingle(sensorRef)
        mi = exposure.getMaskedImage()
        mi /= self.getDarkTime(exposure)
        return exposure

    def getDarkTime(self, exposure):
        if self.config.darkTime is not None:
            return exposure.getMetadata().get(self.config.darkTime)
        return exposure.getCalib().getExpTime()


def flatOverrides(config):
    config.isr.doFlat = False
#   config.isr.doFringe = False


class FlatConfig(DetrendConfig):
    iterations = Field(dtype=int, default=10, doc="Number of iterations for scale determination")
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")

class FlatTask(DetrendTask):
    ConfigClass = FlatConfig
    _DefaultName = "flat"
    overrides = [flatOverrides]
    calibName = "flat"

    def __init__(self, **kwargs):
        super(FlatTask, self).__init__(**kwargs)
        self.makeSubtask("stats")

    def processResult(self, exposure):
        return self.stats.run(exposure)

    def scale(self, ccdKeys, ccdIdLists, data):
        # Format background measurements into a matrix
        indices = dict((name, i) for i, name in enumerate(ccdIdLists.keys()))
        bgMatrix = numpy.array([[0] * len(expList) for expList in ccdIdLists.values()])
        for name in ccdIdLists.keys():
            i = indices[name]
            bgMatrix[i] = data[name]

        self.log.info("Input backgrounds: %s" % bgMatrix)

        # Flat-field scaling
        numCcds = len(ccdIdLists)
        bgMatrix = numpy.log(bgMatrix)      # log(Background) for each exposure/component
        compScales = numpy.zeros(numCcds) # Initial guess at log(scale) for each component
        expScales = numpy.apply_along_axis(lambda x: numpy.average(x - compScales), 0, bgMatrix)

        for iterate in range(self.config.iterations):
            # XXX use masks for each quantity: maskedarrays
            compScales = numpy.apply_along_axis(lambda x: numpy.average(x - expScales), 1, bgMatrix)
            expScales = numpy.apply_along_axis(lambda x: numpy.average(x - compScales), 0, bgMatrix)
            avgScale = numpy.average(numpy.exp(compScales))
            compScales -= numpy.log(avgScale)
            self.log.logdebug("Iteration %d exposure scales: %s" % (iterate, numpy.exp(expScales)))
            self.log.logdebug("Iteration %d component scales: %s" % (iterate, numpy.exp(compScales)))

        expScales = numpy.apply_along_axis(lambda x: numpy.average(x - compScales), 0, bgMatrix)

        if numpy.any(numpy.isnan(expScales)):
            raise RuntimeError("Bad exposure scales: %s --> %s" % (bgMatrix, expScales))

        expScales = numpy.exp(expScales)
        compScales = numpy.exp(compScales)

        self.log.info("Exposure scales: %s" % expScales)
        self.log.info("Component relative scaling: %s" % compScales)

        return dict((ccdName, Struct(ccdScale=compScales[indices[ccdName]], expScales=expScales))
                    for ccdName in ccdIdLists.keys())


class FringeConfig(DetrendConfig):
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")

def fringeOverrides(config):
    pass
#    config.isr.doFringe = False

class FringeTask(DetrendTask):
    """Fringe construction task

    XXX Would like to have this do PCA and generate multiple images, but
    that will take a bit of work with the persistence code.
    """
    ConfigClass = FringeConfig
    _DefaultName = "fringe"
    overrides = [fringeOverrides]
    calibName = "fringe"

    def __init__(self, **kwargs):
        super(FringeTask, self).__init__(**kwargs)
        self.makeSubtask("detection")
        self.makeSubtask("stats")

    def processSingle(self, sensorRef):
        exposure = super(FringeTask, self).processSingle(sensorRef)
        bgLevel = self.stats.run(exposure)
        self.subtractBackground(exposure)
        mi = exposure.getMaskedImage()
        mi /= bgLevel
        footprintSets = self.detection.detectFootprints(exposure)
        mask = exposure.getMaskedImage().getMask()
        detected = mask.addMaskPlane("DETECTED")
        for fpSet in (footprintSets.positive, footprintSets.negative):
            if fpSet is not None:
                afwDet.setMaskFromFootprintList(mask, fpSet.getFootprints(), detected)
        return exposure

    def subtractBackground(self, exposure):
        mi = exposure.getMaskedImage()
        background = measAlg.getBackground(mi, self.config.background).getImageF()
        mi -= background

def getSize(dimList):
    dim = set((w, h) for w,h in dimList)
    dim.discard(None)
    if len(dim) != 1:
        raise RuntimeError("Inconsistent dimensions: %s" % dim)
    return dim.pop()


class MaskCombineConfig(Config):
    maskFraction = Field(doc="Minimum fraction of images where bad pixels got flagged", dtype=float,
                         default=0.5)
    maskPlane = Field(doc="Name of mask plane to set", dtype=str, default="BAD")


class MaskCombineTask(Task):
    ConfigClass = MaskCombineConfig

    def run(footprintSetsList, dimList):
        width, height = getSize(dimList)
        combined = afwImage.MaskU(width, height)
        for footprintSets in footprintSetsList:
            mask = afwImage.MaskU(width, height)
            afwImage.setMaskFromFootprintList(mask, footprintSets.positive, 1)
            afwImage.setMaskFromFootprintList(mask, footprintSets.negative, 1)
            combined += mask

        threshold = afwDet.createThreshold(int(self.config.maskFraction * len(footprintSetsList)))
        footprints = afwDet.FootprintSet(combined, threshold)
        mask = combined.addMaskPlane(self.config.maskPlane)
        combined.set(0)
        combined.setMaskFromFootprintList(footprints, mask)

        return combined


class MaskConfig(DetrendConfig):
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    def setDefaults(self):
        self.combination.retarget(MaskCombineTask)

def maskOverrides(config):
    pass

class MaskTask(DetrendTask):
    ConfigClass = MaskConfig
    _DefaultName = "mask"
    overrides = [maskOverrides]
    calibName = "mask"

    def __init__(self, **kwargs):
        super(MaskTask, self).__init__(**kwargs)
        background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
        self.makeSubtask("detection")

    def processWrite(self, *args, **kwargs):
        """Don't write anything --- there's no need"""
        pass

    def processResult(self, exposure):
        self.subtractBackground(exposure)
        footprintSets = self.detection.detectFootprints(exposure)
        return Struct(dim=exposure.getDimensions(), footprints=footprintSets)

    def subtractBackground(self, exposure):
        mi = exposure.getMaskedImage()
        background = measAlg.getBackground(mi, self.config.background).getImageF()
        mi -= background

    def scales(self, data):
        return data

    def combine(self, struct):
        fpSetsList = [s.footprints for s in struct.scales]
        dimsList = [s.dim for s in struct.scales]
        combined = self.combination.run(fpSetsList, dimsList)
        self.write(combined, struct.outputId)
