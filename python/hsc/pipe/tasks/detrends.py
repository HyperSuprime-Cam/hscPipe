#!/usr/bin/env python

import os
import sys
import numpy
import argparse
import traceback

from lsst.pex.config import Config, ConfigField, ConfigurableField, Field, ListField
from lsst.pipe.base import Task, Struct, TaskRunner
import lsst.daf.base as dafBase
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom
import lsst.meas.algorithms as measAlg
import lsst.afw.geom.ellipses as afwEll

import lsst.obs.subaru.isr as hscIsr

import hsc.pipe.base.butler as hscButler
from hsc.pipe.base.mpi import MpiTask, MpiArgumentParser, abortOnError, thisNode
from hsc.pipe.base.pbs import PbsCmdLineTask

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
    """Configuration for combining detrend images"""
    rows = Field(doc="Number of rows to read at a time", dtype=int, default=128)
    mask = ListField(doc="Mask planes to respect", dtype=str, default=["SAT", "DETECTED", "INTRP"])
    combine = Field(doc="Statistic to use for combination (from lsst.afw.math)", dtype=int,
                    default=afwMath.MEANCLIP)
    clip = Field(doc="Clipping threshold for combination", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for combination", dtype=int, default=3)
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")

class DetrendCombineTask(Task):
    """Task to combine detrend images"""
    ConfigClass = DetrendCombineConfig

    def __init__(self, *args, **kwargs):
        super(DetrendCombineTask, self).__init__(*args, **kwargs)
        self.makeSubtask("stats")

    def run(self, sensorRefList, expScales=None, finalScale=None, inputName="postISRCCD"):
        """Combine detrend images for a single sensor

        @param sensorRefList   List of data references to combine (for a single sensor)
        @param expScales       List of scales to apply for each exposure
        @param finalScale      Desired scale for final combined image
        @param inputName       Data set name for inputs
        @return combined image
        """
        width, height = self.getDimensions(sensorRefList)
        maskVal = 0
        for mask in self.config.mask:
            maskVal |= afwImage.MaskU.getPlaneBitMask(mask)
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
                exposure = sensorRef.get(inputName + "_sub", bbox=box)
                if expScales is not None:
                    self.applyScale(exposure, expScales[i])
                imageList[i] = exposure.getMaskedImage()

            self.combine(subCombined, imageList, stats)

        if finalScale is not None:
            background = self.stats.run(combined)
            self.log.info("%s: Measured background of stack is %f; adjusting to %f" %
                         (thisNode(), background, finalScale))
            combined *= finalScale / background

        return combined

    def getDimensions(self, sensorRefList, inputName="postISRCCD"):
        """Get dimensions of the inputs"""
        dimList = []
        for sensorRef in sensorRefList:
            md = sensorRef.get(inputName + "_md")
            dimList.append(afwGeom.Extent2I(md.get("NAXIS1"), md.get("NAXIS2")))
        return getSize(dimList)

    def applyScale(self, exposure, scale=None):
        """Apply scale to input exposure"""
        if scale is not None:
            mi = exposure.getMaskedImage()
            mi /= scale

    def combine(self, target, imageList, stats):
        """Combine multiple images

        @param target      Target image to receive the combined pixels
        @param imageList   List of input images
        @param stats       Statistics control
        """
        if False:
            # In-place stacks are now supported on LSST's afw, but not yet on HSC
            afwMath.statisticsStack(target, imageList, self.config.combine, stats)
        else:
            stack = afwMath.statisticsStack(imageList, self.config.combine, stats)
            target <<= stack.getImage()



def getCcdName(ccdId, ccdKeys):
    """Return the 'CCD name' from the data identifier

    The 'CCD name' is a tuple of the values in the data identifier
    that identify the CCD.  The names in the data identifier that
    identify the CCD is provided as 'ccdKeys'.

    @param ccdId    Data identifier for CCD
    @param ccdKeys  Data identifier keys for the 'sensor' level
    @return ccd name
    """
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

    @param expRefList   List of data references for exposures
    @param level        Level for the butler to generate CCDs
    @return CCD keywords, dict of data identifier lists for each CCD
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
    """Split name=value pairs and put the result in a dict"""
    def __call__(self, parser, namespace, values, option_string):
        output = getattr(namespace, self.dest, {})
        for nameValue in values:
            name, sep, valueStr = nameValue.partition("=")
            if not valueStr:
                parser.error("%s value %s must be in form name=value" % (option_string, nameValue))
            output[name] = valueStr
        setattr(namespace, self.dest, output)

class DetrendArgumentParser(MpiArgumentParser):
    """Add a --detrendId argument to the argument parser"""
    def __init__(self, calibName, *args, **kwargs):
        super(DetrendArgumentParser, self).__init__(*args, **kwargs)
        self.calibName = calibName
        self.add_id_argument("--id", datasetType="raw", level="visit",
                             help="input identifiers, e.g., --id visit=123 ccd=4", rootOnly=False)
        self.add_argument("--detrendId", nargs="*", action=DetrendIdAction, default={},
                          help="identifiers for detrend, e.g., --detrendId version=1",
                          metavar="KEY=VALUE1[^VALUE2[^VALUE3...]")
    def parse_args(self, *args, **kwargs):
        namespace = super(DetrendArgumentParser, self).parse_args(*args, **kwargs)

        keys = namespace.butler.getKeys(self.calibName)
        parsed = {}
        for name, value in namespace.detrendId.items():
            if not name in keys:
                self.error("%s is not a relevant detrend identifier key (%s)" % (name, keys))
            parsed[name] = keys[name](value)
        namespace.detrendId = parsed

        return namespace

class DetrendConfig(Config):
    """Configuration for constructing detrends"""
    isr = ConfigurableField(target=hscIsr.SubaruIsrTask, doc="ISR configuration")
    dateObs = Field(dtype=str, default="dateObs", doc="Key for observation date in exposure registry")
    dateCalib = Field(dtype=str, default="calibDate", doc="Key for detrend date in calib registry")
    filter = Field(dtype=str, default="filter", doc="Key for filter name in exposure/calib registries")
    combination = ConfigurableField(target=DetrendCombineTask, doc="Detrend combination configuration")
    def setDefaults(self):
        self.isr.doWrite = False

class DetrendTaskRunner(TaskRunner):
    """Get parsed values into the DetrendTask.run"""
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [dict(expRefList=parsedCmd.id.refList, butler=parsedCmd.butler, detrendId=parsedCmd.detrendId)]

    def __call__(self, args):
        task = self.TaskClass(config=self.config, log=self.log)
        if self.doRaise:
            result = task.run(**args)
        else:
            try:
                result = task.run(**args)
            except Exception, e:
                task.log.fatal("Failed: %s" % e)
                traceback.print_exc(file=sys.stderr)

        if self.doReturnResults:
            return Struct(
                args = args,
                metadata = task.metadata,
                result = result,
            )

class DetrendTask(PbsCmdLineTask, MpiTask):
    """Base class for constructing detrends.

    This should be subclassed for each of the required detrend types.
    The subclass should be sure to define the following class variables:
    * _DefaultName: default name of the task, used by CmdLineTask
    * calibName: name of the calibration data set in the butler
    * overrides: a list of functions for setting a configuration, used by CmdLineTask
    """
    ConfigClass = DetrendConfig
    RunnerClass = DetrendTaskRunner

    def __init__(self, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(DetrendTask, self).__init__(**kwargs)
        self.makeSubtask("isr")
        self.makeSubtask("combination")

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numCcds = sum(1 for raft in parsedCmd.butler.get("camera") for ccd in cameraGeom.cast_Raft(raft))
        numExps = len(cls.RunnerClass.getTargetList(parsedCmd)[0]['expRefList'])
        numCycles = int(numCcds/float(numNodes*numProcs) + 0.5)
        return time*numExps*numCycles

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doPbs = kwargs.pop("doPbs", False)
        return DetrendArgumentParser(calibName=cls.calibName, name=cls._DefaultName, *args, **kwargs)

    @abortOnError
    def run(self, expRefList, butler, detrendId):
        """Construct a detrend from a list of exposure references

        This is the entry point, called by the TaskRunner.__call__

        All nodes execute this method.

        @param expRefList  List of data references at the exposure level
        @param butler      Data butler
        @param detrendId   Identifier dict for detrend
        """
        self.butler = butler

        if self.rank == self.root:
            outputId = self.getOutputId(expRefList, detrendId)
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

        import pbasf2
        ccdIdLists = pbasf2.Broadcast(self.comm, ccdIdLists, root=self.root)

        # Scatter: process CCDs independently
        data = self.scatterProcess(ccdKeys, ccdIdLists)

        # Gather: determine scalings
        if self.rank == self.root:
            scales = self.scale(ccdKeys, ccdIdLists, data)
        else:
            scales = None, None

        # Scatter: combine
        self.scatterCombine(outputId, ccdKeys, ccdIdLists, scales)

    def getOutputId(self, expRefList, detrendId):
        """Generate the data identifier for the output detrend

        The mean date and the common filter are included, using keywords
        from the configuration.  The CCD-specific part is not included
        in the data identifier.

        Only the root node executes this method (it will share the results with the slaves).

        @param expRefList  List of data references at exposure level
        """
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
        outputId.update(detrendId)
        return outputId

    def getMjd(self, dataId):
        """Determine the Modified Julian Date (MJD) from a data identifier"""
        dateObs = dataId[self.config.dateObs]
        try:
            dt = dafBase.DateTime(dateObs)
        except:
            dt = dafBase.DateTime(dateObs + "T12:00:00.0Z")
        return dt.get(dafBase.DateTime.MJD)

    def getFilter(self, dataId):
        """Determine the filter from a data identifier"""
        return dataId[self.config.filter]

    def scatterProcess(self, ccdKeys, ccdIdLists):
        """Scatter the processing among the nodes

        We scatter the data wider than the just the number of CCDs, to make
        full use of all available processors.  This necessitates piecing
        everything back together in the same format as ccdIdLists afterwards.

        All nodes execute this method (though with different behaviour).

        @param ccdKeys     Keywords that identify a CCD in the data identifier
        @param ccdIdLists  Dict of data identifier lists for each CCD name
        @return Dict of lists of returned data for each CCD name
        """
        if self.rank == self.root:
            dataIdList = sum(ccdIdLists.values(), [])
        else:
            dataIdList = None

        self.log.info("Scatter processing on %s" % thisNode())
        import pbasf2
        resultList = pbasf2.ScatterJob(self.comm, self.process, dataIdList, root=self.root)
        if self.rank == self.root:
            # Piece everything back together
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
        """Process a CCD, specified by a data identifier

        Only slave nodes execute this method.
        """
        self.log.info("Processing %s on %s" % (ccdId, thisNode()))
        sensorRef = hscButler.getDataRef(self.butler, ccdId)
        exposure = self.processSingle(sensorRef)
        self.processWrite(sensorRef, exposure)
        return self.processResult(exposure)

    def processSingle(self, dataRef):
        """Process a single CCD, specified by a data reference

        Only slave nodes execute this method.
        """
        return self.isr.run(dataRef).exposure

    def processWrite(self, dataRef, exposure, outputName="postISRCCD"):
        """Write the processed CCD

        Only slave nodes execute this method.

        @param dataRef     Data reference
        @param exposure    CCD exposure to write
        @param outputName  Data type name for butler.
        """
        dataRef.put(exposure, outputName)

    def processResult(self, exposure):
        """Extract processing results from a processed exposure

        Only slave nodes execute this method.  This method generates
        what is gathered by the master node --- it must be picklable!
        """
        return None

    def scale(self, ccdKeys, ccdIdLists, data):
        """Determine scaling across CCDs and exposures

        This is necessary mainly for flats, so as to determine a
        consistent scaling across the entire focal plane.

        Only the master node executes this method.

        @param ccdKeys     Keywords that identify a CCD in the data identifier
        @param ccdIdLists  Dict of data identifier lists for each CCD name
        @param data        Dict of lists of returned data for each CCD name
        @return dict of Struct(ccdScale: scaling for CCD,
                               expScales: scaling for each exposure
                               ) for each CCD name
        """
        self.log.info("Scale on %s" % thisNode())
        return dict((name, Struct(ccdScale=None, expScales=[None] * len(ccdIdLists[name])))
                    for name in ccdIdLists.keys())

    def scatterCombine(self, outputId, ccdKeys, ccdIdLists, scales):
        """Scatter the combination across multiple nodes

        In this case, we can only scatter across as many nodes as
        there are CCDs.

        All nodes execute this method (though with different behaviour).

        @param outputId    Output identifier (exposure part only)
        @param ccdKeys     Keywords that identify a CCD in the data identifier
        @param ccdIdLists  Dict of data identifier lists for each CCD name
        @param scales      Dict of structs with scales, for each CCD name
        """
        self.log.info("Scatter combination on %s" % thisNode())
        if self.rank == self.root:
            data = [Struct(ccdIdList=ccdIdLists[ccdName], scales=scales[ccdName],
                           outputId=dict(outputId.items() + [(k,ccdName[i]) for i, k in enumerate(ccdKeys)]))
                    for ccdName in ccdIdLists.keys()]
        else:
            data = None
        import pbasf2
        pbasf2.ScatterJob(self.comm, self.combine, data, root=self.root)

    def combine(self, struct):
        """Combine multiple exposures of a particular CCD and write the output

        Only the slave nodes execute this method.

        The input is a struct consisting of the following components:
        @param ccdIdList   List of data identifiers for combination
        @param scales      Scales to apply (expScales are scalings for each exposure,
                           ccdScale is final scale for combined image)
        @param outputId    Data identifier for combined image (fully qualified for this CCD)
        """
        dataRefList = [hscButler.getDataRef(self.butler, dataId) for dataId in struct.ccdIdList]
        self.log.info("Combining %s on %s" % (struct.outputId, thisNode()))
        detrend = self.combination.run(dataRefList, expScales=struct.scales.expScales,
                                       finalScale=struct.scales.ccdScale)
        self.write(detrend, struct.outputId)

    def write(self, exposure, dataId):
        """Write the final combined detrend

        Only the slave nodes execute this method

        @param exposure  CCD exposure to write
        @param dataId    Data identifier
        """
        self.log.info("Writing %s on %s" % (dataId, thisNode()))
        self.butler.put(exposure, self.calibName, dataId)


class BiasConfig(DetrendConfig):
    """Configuration for bias construction.

    No changes required compared to the base class, but
    subclassed for distinction.
    """
    pass


class BiasTask(DetrendTask):
    """Bias construction"""
    ConfigClass = BiasConfig
    _DefaultName = "bias"
    calibName = "bias"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for bias construction"""
        config.isr.doBias = False
        config.isr.doDark = False
        config.isr.doFlat = False
#        config.isr.doFringe = False


class DarkConfig(DetrendConfig):
    """Configuration for dark construction"""
    darkTime = Field(dtype=str, default="DARKTIME", doc="Header keyword for time since last CCD wipe")

class DarkTask(DetrendTask):
    """Dark construction

    The only major difference from the base class is dividing
    each image by the dark time to generate images of the
    dark rate.
    """
    ConfigClass = DarkConfig
    _DefaultName = "dark"
    calibName = "dark"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides to apply for dark construction"""
        config.isr.doDark = False
        config.isr.doFlat = False
#        config.isr.doFringe = False

    def processSingle(self, sensorRef):
        """Divide each processed image by the dark time to generate images of the dark rate"""
        exposure = super(DarkTask, self).processSingle(sensorRef)
        mi = exposure.getMaskedImage()
        mi /= self.getDarkTime(exposure)
        return exposure

    def getDarkTime(self, exposure):
        """Retrieve the dark time"""
        if self.config.darkTime is not None:
            return exposure.getMetadata().get(self.config.darkTime)
        return exposure.getCalib().getExpTime()


class FlatCombineConfig(DetrendCombineConfig):
    """Configuration for flat construction"""
    doJacobian = Field(dtype=bool, default=False, doc="Apply Jacobian to flat-field?")


class FlatCombineTask(DetrendCombineTask):
    """Combination for flat-fields

    We allow the flat-field to be corrected for the Jacobian.

    The observed flat-field has a constant exposure per unit area.
    However, distortion in the camera makes the pixel area (the angle
    on the sky subtended per pixel) larger as one moves away from the
    optical axis, so that the flux in the observed flat-field drops.
    But this drop does not mean the detector is less sensitive.  The
    Jacobian is a rough estimate of the relative area of each pixel.
    The correction is therefore achieved by multiplying the observed
    flat-field by the Jacobian to create a "photometric flat" which
    has constant exposure per pixel --- a true measure of the
    point-source sensitivity of the camera as a function of pixel
    (modulo contributions from scattered light, which require much
    more work to account for).

    Note, however, that application of this correction means that
    images flattened with this (Jacobian-corrected) "photometric flat"
    will not have a flat sky, potentially making sky subtraction more
    difficult.  Furthermore, care must be taken to ensure this
    correction is not applied more than once (e.g., in warping).
    """
    ConfigClass = FlatCombineConfig

    def run(self, sensorRefList, expScales=None, finalScale=None):
        """Multiply the combined flat-field by the Jacobian"""
        combined = super(FlatCombineTask, self).run(sensorRefList, expScales=expScales, finalScale=finalScale)
        if self.config.doJacobian:
            jacobian = self.getJacobian(sensorRefList[0], combined.getDimensions())
            combined *= jacobian
        return combined

    def getJacobian(self, sensorRef, dimensions, inputName="postISRCCD"):
        """Calculate the Jacobian as a function of position

        @param sensorRef    Data reference for a representative CCD (to get the Detector)
        @param dimensions   Dimensions of the flat-field
        @param inputName    Data set name for inputs
        @return Jacobian image
        """
        # Retrieve the detector and distortion
        # XXX It's unfortunate that we have to read an entire image to get the detector, but there's no
        # public API in the butler to get the same.
        image = sensorRef.get(inputName)
        detector = image.getDetector()
        distortion = detector.getDistortion()
        del image

        # Calculate the Jacobian for each pixel
        # XXX This would be faster in C++, but it's not awfully slow.
        jacobian = afwImage.ImageF(dimensions)
        array = jacobian.getArray()
        width, height = dimensions
        circle = afwEll.Quadrupole(1.0, 1.0, 0.0)
        for y in xrange(height):
            for x in xrange(width):
                array[y,x] = distortion.distort(afwGeom.Point2D(x, y), circle, detector).getDeterminant()

        return jacobian


class FlatConfig(DetrendConfig):
    """Configuration for flat construction"""
    iterations = Field(dtype=int, default=10, doc="Number of iterations for scale determination")
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")
    def setDefaults(self):
        self.combination.retarget(FlatCombineTask)

class FlatTask(DetrendTask):
    """Flat construction

    The principal change involves gathering the background values from each
    image and using them to determine the scalings for the final combination.
    """
    ConfigClass = FlatConfig
    _DefaultName = "flat"
    calibName = "flat"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides for flat construction"""
        config.isr.doFlat = False
#        config.isr.doFringe = False


    def __init__(self, **kwargs):
        super(FlatTask, self).__init__(**kwargs)
        self.makeSubtask("stats")

    def processResult(self, exposure):
        return self.stats.run(exposure)

    def scale(self, ccdKeys, ccdIdLists, data):
        """Determine the scalings for the final combination

        We have a matrix B_ij = C_i E_j, where C_i is the relative scaling
        of one CCD to all the others in an exposure, and E_j is the scaling
        of the exposure.  We determine the C_i and E_j from B_ij by iteration,
        under the additional constraint that the average CCD scale is unity.
        We convert everything to logarithms so we can work linearly.  This
        algorithm comes from Eugene Magnier and Pan-STARRS.
        """
        # Format background measurements into a matrix
        indices = dict((name, i) for i, name in enumerate(ccdIdLists.keys()))
        bgMatrix = numpy.array([[0] * len(expList) for expList in ccdIdLists.values()])
        for name in ccdIdLists.keys():
            i = indices[name]
            bgMatrix[i] = data[name]

        numpyPrint = numpy.get_printoptions()
        numpy.set_printoptions(threshold='nan')
        self.log.info("Input backgrounds: %s" % bgMatrix)

        # Flat-field scaling
        numCcds = len(ccdIdLists)
        numExps = bgMatrix.shape[1]
        bgMatrix = numpy.log(bgMatrix)      # log(Background) for each exposure/component
        bgMatrix = numpy.ma.masked_array(bgMatrix, numpy.isnan(bgMatrix))
        compScales = numpy.zeros(numCcds) # Initial guess at log(scale) for each component
        expScales = numpy.array([(bgMatrix[:,i] - compScales).mean() for i in range(numExps)])

        for iterate in range(self.config.iterations):
            compScales = numpy.array([(bgMatrix[i,:] - expScales).mean() for i in range(numCcds)])
            expScales = numpy.array([(bgMatrix[:,i] - compScales).mean() for i in range(numExps)])

            avgScale = numpy.average(numpy.exp(compScales))
            compScales -= numpy.log(avgScale)
            self.log.logdebug("Iteration %d exposure scales: %s" % (iterate, numpy.exp(expScales)))
            self.log.logdebug("Iteration %d component scales: %s" % (iterate, numpy.exp(compScales)))

        expScales = numpy.array([(bgMatrix[:,i] - compScales).mean() for i in range(numExps)])

        if numpy.any(numpy.isnan(expScales)):
            raise RuntimeError("Bad exposure scales: %s --> %s" % (bgMatrix, expScales))

        expScales = numpy.exp(expScales)
        compScales = numpy.exp(compScales)

        self.log.info("Exposure scales: %s" % expScales)
        self.log.info("Component relative scaling: %s" % compScales)

        return dict((ccdName, Struct(ccdScale=compScales[indices[ccdName]], expScales=expScales))
                    for ccdName in ccdIdLists.keys())


class FringeConfig(DetrendConfig):
    """Configuration for fringe construction"""
    stats = ConfigurableField(target=DetrendStatsTask, doc="Background statistics configuration")
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    detectSigma = Field(dtype=float, default=1.0, doc="Detection PSF gaussian sigma")

class FringeTask(DetrendTask):
    """Fringe construction task

    The principal change from the base class is that the images are
    background-subtracted and rescaled by the background.

    XXX This is probably not right for a straight-up combination, as we
    are currently doing, since the fringe amplitudes need not scale with
    the continuum.

    XXX Would like to have this do PCA and generate multiple images, but
    that will take a bit of work with the persistence code.
    """
    ConfigClass = FringeConfig
    _DefaultName = "fringe"
    calibName = "fringe"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides for fringe construction"""
        config.isr.doFringe = False

    def __init__(self, **kwargs):
        super(FringeTask, self).__init__(**kwargs)
        self.makeSubtask("detection")
        self.makeSubtask("stats")

    def processSingle(self, sensorRef):
        """Subtract the background and normalise by the background level"""
        exposure = super(FringeTask, self).processSingle(sensorRef)
        bgLevel = self.stats.run(exposure)
        self.subtractBackground(exposure)
        mi = exposure.getMaskedImage()
        mi /= bgLevel
        footprintSets = self.detection.detectFootprints(exposure, sigma=self.config.detectSigma)
        mask = exposure.getMaskedImage().getMask()
        detected = mask.addMaskPlane("DETECTED")
        for fpSet in (footprintSets.positive, footprintSets.negative):
            if fpSet is not None:
                afwDet.setMaskFromFootprintList(mask, fpSet.getFootprints(), detected)
        return exposure

    def subtractBackground(self, exposure):
        """Subtract the background from the provided exposure"""
        mi = exposure.getMaskedImage()
        background = measAlg.getBackground(mi, self.config.background).getImageF()
        mi -= background

def getSize(dimList):
    """Determine the consistent size, given a list of image sizes"""
    dim = set((w, h) for w,h in dimList)
    dim.discard(None)
    if len(dim) != 1:
        raise RuntimeError("Inconsistent dimensions: %s" % dim)
    return dim.pop()


class MaskCombineConfig(Config):
    """Configuration for generating a mask from individual footprints of discrepant pixels"""
    maskFraction = Field(doc="Minimum fraction of images where bad pixels got flagged", dtype=float,
                         default=0.5)
    maskPlane = Field(doc="Name of mask plane to set", dtype=str, default="BAD")


class MaskCombineTask(Task):
    """Generate a mask from individual footprints of discrepant pixels"""
    ConfigClass = MaskCombineConfig

    def run(footprintSetsList, dimList):
        """Generate a mask

        We are given a set of footprints of discrepant pixels for each image,
        and identify bad pixels as those which are consistently discrepant.

        @param footprintSetsList    List of footprintSets for each image
        @param dimList              List of sizes for each image
        @return Combined image
        """
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
    """Configuration for mask construction"""
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    def setDefaults(self):
        self.combination.retarget(MaskCombineTask)

class MaskTask(DetrendTask):
    """Mask construction task"""
    ConfigClass = MaskConfig
    _DefaultName = "mask"
    calibName = "mask"

    @classmethod
    def applyOverrides(cls, config):
        """Overrides for mask construction"""
        pass

    def __init__(self, **kwargs):
        super(MaskTask, self).__init__(**kwargs)
        background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")
        self.makeSubtask("detection")

    def processWrite(self, *args, **kwargs):
        """Don't write anything

        There's no need to write anything, as all the required information
        is returned.
        """
        pass

    def processResult(self, exposure):
        """Return the discrepant pixels

        XXX Not sure if footprintSets will pickle
        """
        self.subtractBackground(exposure)
        footprintSets = self.detection.detectFootprints(exposure)
        return Struct(dim=exposure.getDimensions(), footprints=footprintSets)

    def subtractBackground(self, exposure):
        """Subtract background from exposure"""
        mi = exposure.getMaskedImage()
        background = measAlg.getBackground(mi, self.config.background).getImageF()
        mi -= background

    def scales(self, data):
        """No scaling required.

        We just want to propagate the data returned from processing.
        """
        return data

    def combine(self, struct):
        """Combine the lists of discrepant pixels"""
        fpSetsList = [s.footprints for s in struct.scales]
        dimsList = [s.dim for s in struct.scales]
        combined = self.combination.run(fpSetsList, dimsList)
        self.write(combined, struct.outputId)
