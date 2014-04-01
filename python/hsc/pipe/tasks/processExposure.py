import sys
import math
import itertools
import collections

import hsc.pipe.tasks.plotSetup
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCg
import hsc.pipe.base.butler as hscButler
from lsst.pipe.base import Struct, ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField, ConfigField
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask
from hsc.pipe.tasks.photometricSolution import PhotometricSolutionTask
from hsc.pipe.tasks.backgroundModels import WarpedBackground, BackgroundConfig
from hsc.pipe.base.pool import abortOnError, NODE, Pool, Debugger
from hsc.pipe.base.pbs import PbsPoolTask
from hsc.meas.tansip.solvetansip import SolveTansipTask

Debugger().enabled = True


class ProcessExposureConfig(Config):
    processCcd = ConfigurableField(target=SubaruProcessCcdTask, doc="CCD processing task")
    photometricSolution = ConfigurableField(target=PhotometricSolutionTask, doc="Global photometric solution")
    instrument = Field(dtype=str, default="suprimecam", doc="Instrument name, for solvetansip")
    doSolveTansip = Field(dtype=bool, default=True, doc="Run solvetansip?")
    solveTansip = ConfigurableField(target=SolveTansipTask, doc="Global astrometric solution")
    doPhotometricSolution = Field(dtype=bool, default=False, doc="Run global photometric solution?")
    doBackground = Field(dtype=bool, default=False, doc="Do global background subtraction?")
    background = ConfigField(dtype=BackgroundConfig, doc="Global background subtraction config")

    def setDefaults(self):
        # We will do persistence ourselves
        self.processCcd.isr.doWrite = False
        self.processCcd.doWriteCalibrate = False
        self.processCcd.doWriteSources = False
        self.processCcd.doWriteHeavyFootprintsInSources = False
        self.processCcd.doFinalWrite = False

        self.background.xSize = 2048
        self.background.ySize = 2048
        self.background.extrapolationMode = "FULL"


class ProcessExposureTask(PbsPoolTask):
    """Process an entire exposure at once.

    We use MPI to gather the match lists for exposure-wide astrometric and
    photometric solutions.  Note that because of this, different nodes
    see different parts of the code.
    """

    RunnerClass = hscButler.ButlerTaskRunner
    ConfigClass = ProcessExposureConfig
    _DefaultName = "processExposure"

    def __init__(self, *args, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessExposureTask, self).__init__(*args, **kwargs)
        self.makeSubtask("processCcd")
        self.makeSubtask("photometricSolution", schema=self.processCcd.schema)
        self.makeSubtask("solveTansip")

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numCcds = sum(1 for raft in parsedCmd.butler.get("camera") for ccd in afwCg.cast_Raft(raft))
        numCycles = int(math.ceil(numCcds/float(numNodes*numProcs)))
        numExps = len(cls.RunnerClass.getTargetList(parsedCmd))
        return time*numExps*numCycles

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doPbs = kwargs.pop("doPbs", False)
        parser = ArgumentParser(name="processExposure", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", level="visit",
                               help="data ID, e.g. --id visit=12345")
        return parser

    @abortOnError
    def run(self, expRef, butler):
        """Process a single exposure, with scatter-gather-scatter using MPI.

        All nodes execute this method, though the master and slaves have different
        routes through it.  The expRef is only a DummyDataRef on the slaves.
        """
        pool = Pool("processExposure")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        dataIdList = dict([(ccdRef.get("ccdExposureId"), ccdRef.dataId)
                           for ccdRef in expRef.subItems("ccd") if ccdRef.datasetExists("raw")])
        dataIdList = collections.OrderedDict(sorted(dataIdList.items()))

        # Scatter: process CCDs independently
        structList = pool.map(self.process, dataIdList.values())
        numGood = sum(1 for s in structList if s is not None)
        if numGood == 0:
            self.log.warn("All CCDs in exposure failed")
            return

        # Gathered: global WCS solution
        matchLists = self.getMatchLists(dataIdList, structList)
        wcsList = self.solveAstrometry(matchLists, butler.mapper.camera)
        fluxMag0 = self.solvePhotometry(matchLists.values(), self.getFilterName(structList))
        ccdIdList = dataIdList.keys()
        if self.config.doBackground:
            bg = self.backgroundPrepare(wcsList, structList)
            bgList = pool.mapToPrevious(self.measureBackground, dataIdList.values(), bg)
            bg = self.backgroundReduce(bg, bgList)
        else:
            bg = None

        # Scatter with data from root: save CCDs with update astrometric/photometric solutions
        solutionList = [Struct(ccdId=ccdId, fluxMag0=fluxMag0, dataId=dataIdList[ccdId],
                               wcs=wcsList[ccdId] if ccdId in wcsList else None)
                        for ccdId in ccdIdList]
        pool.mapToPrevious(self.write, solutionList, bg)

    def process(self, cache, dataId):
        """Process a single CCD and save the results for a later write.

        Only slaves execute this method.
        """
        cache.result = None
        dataRef = hscButler.getDataRef(cache.butler, dataId)
        ccdId = dataRef.get("ccdExposureId")
        self.log.info("Started processing %s (ccdId=%d) on %s" % (dataId, ccdId, NODE))
        try:
            result = self.processCcd.run(dataRef)
        except Exception, e:
            self.log.warn("Failed to process %s: %s\n" % (dataId, e))
            import traceback
            traceback.print_exc()
            return None

        # Cache the results (in particular, the image)
        cache.result = result
        filterName = result.exposure.getFilter().getName()

        # Reformat the matches for MPI transfer
        matches, numMatches = None, 0
        box = None
        if result.matches is not None and result.matchMeta is not None:
            matches = afwTable.ReferenceMatchVector(result.matches)
            numMatches = len(matches)
            box = result.exposure.getMaskedImage().getBBox()

        self.log.info("Finished processing %s (ccdId=%d) on %s with %d matches" %
                      (dataId, ccdId, NODE, numMatches))

        return Struct(ccdId=ccdId, matches=matches, filterName=filterName, box=box)

    def getFilterName(self, structList):
        """Determine the filter name from the list of structs returned by process().

        Only the master executes this method, as the structList is only valid there.
        """
        filterList = [s.filterName if s is not None else None for s in structList]
        filterSet = set(filterList)
        filterSet.discard(None) # Just in case
        if len(filterSet) != 1:
            raise RuntimeError("Multiple filters over exposure: %s" % filterSet)
        return filterSet.pop()

    def getMatchLists(self, dataIdList, structList):
        """Generate a list of matches for each CCD from the list of structs returned by process().

        The matches are reconsituted from the transfer format.

        Only the master executes this method, as the structList is only valid there.
        """
        lookup = dict((s.ccdId, s.matches) for s in structList if s is not None)
        return collections.OrderedDict((ccdId, lookup.get(ccdId, None)) for ccdId in dataIdList)

    def solveAstrometry(self, matchLists, cameraGeom):
        """Determine a global astrometric solution for the exposure.

        Only the master executes this method, as the matchLists is only valid there.
        """

        wcsList = [None] * len(matchLists)
        if self.config.doSolveTansip:
            try:
                solvetansipIn = [self.solveTansip.convert(ml) if ml is not None else []
                                 for ml in matchLists.itervalues()]

                wcsList = self.solveTansip.solve(self.config.instrument, cameraGeom, solvetansipIn)
            except Exception, e:
                self.log.warn("WARNING: Global astrometric solution failed: %s\n" % e)
        else:
            self.log.info("solvetansip disabled in configuration")
        return dict(zip(matchLists.keys(), wcsList))

    def solvePhotometry(self, matchLists, filterName):
        if not self.config.doPhotometricSolution:
            return None
        try:
            return self.photometricSolution.run(matchLists, filterName)
        except Exception, e:
            self.log.warn("Failed to determine global photometric zero-point: %s" % e)
        return None

    def backgroundPrepare(self, wcsDict, structList):
        """Prepare to measure global background

        We generate a consistent coordinate system over the exposure,
        and use that to seed a background model which can be used to
        make measurements on the individual CCDs.

        This method is run only on the master node.

        @param wcsDict: Dict of ccdId --> WCS solutions
        @param structList: List of process results
        @return virgin background model
        """
        # Determine bounds of images
        xMin, xMax = 2, -2
        yMin, yMax = 2, -2
        zMin, zMax = 2, -2
        scale = 0.0
        num = 0
        for struct in structList:
            if struct is None:
                continue
            box = struct.box
            if box is None:
                continue
            ccdId = struct.ccdId
            if not ccdId in wcsDict:
                continue
            wcs = wcsDict[struct.ccdId]
            if wcs is None:
                continue
            scale += wcs.pixelScale().asDegrees()
            for p in box.getCorners():
                xyz = wcs.pixelToSky(afwGeom.Point2D(p)).getVector()
                xMin = min(xMin, xyz.getX())
                xMax = max(xMax, xyz.getX())
                yMin = min(yMin, xyz.getY())
                yMax = max(yMax, xyz.getY())
                zMin = min(zMin, xyz.getZ())
                zMax = max(zMax, xyz.getZ())

            num += 1
        # Create Wcs
        center = afwGeom.Point3D(0.5*(xMin+xMax), 0.5*(yMin+yMax), 0.5*(zMin+xMax))
        scale /= num
        wcs = afwImage.makeWcs(afwCoord.Coord(xyz), afwGeom.Point2D(0, 0), scale, 0.0, 0.0, scale)
        # Convert sky xyz bounds to image xy box
        box = afwGeom.Box2I()
        for x, y, z in itertools.product((xMin, xMax), (yMin, yMax), (zMin, zMax)):
            box.include(afwGeom.Point2I(wcs.skyToPixel(afwCoord.Coord(afwGeom.Point3D(x, y, z)))))

        return WarpedBackground(self.config.background, box, wcs)

    def measureBackground(self, cache, dataId, bg):
        """Measure background on individual CCD

        The global background model is used to measure on the CCD.
        The resultant measurement is returned in order to form a
        populated full-exposure background model.

        Only slave nodes run this method.

        @param cache: Pool cache
        @param dataId: Data identifier; unused (only required for Pool API)
        @param bg: Unpopulated global background model
        @return Global background model populated with CCD measurements
        """
        if not cache.result:
            return None
        # Reinstate old background measurement
        exposure = cache.result.exposure
        if not exposure:
            return None
        exposure.getMaskedImage().__iadd__(cache.result.backgrounds.getImage())
        # Measure new background
        bg.addImage(exposure)
        return bg

    def backgroundReduce(self, bgMaster, bgList):
        """Consolidate background models

        Reduce the individual CCD measurements of the global
        background model to produce a full-exposure background
        model.

        @param bgMaster: Master (unpopulated) background model
        @param bgList: List of background models from each CCD
        @return Fully-populated background model
        """
        bgList = [bg for bg in bgList if bg is not None]
        return reduce(lambda x,y: x.merge(y), bgList, bgMaster)

    def write(self, cache, struct, bg):
        """Write the outputs.

        The cached results are written along with revised astrometric and photometric solutions.

        This method is only executed on the slaves.
        """
        if not cache.result or not struct:
            # Processing must have failed: nothing we can do
            return

        if bg is not None:
            image = cache.result.exposure.getMaskedImage()
            image -= bg.getImage(image.getBBox(), struct.wcs)

        dataId = struct.dataId
        ccdId = struct.ccdId
        dataRef = hscButler.getDataRef(cache.butler, dataId)
        self.log.info("Start writing %s (ccdId=%d) on %s" % (dataId, ccdId, NODE))
        wcs = struct.wcs
        fluxMag0 = struct.fluxMag0

        try:
            self.processCcd.write(dataRef, cache.result, wcs=wcs, fluxMag0=fluxMag0)
            del cache.result
        except Exception, e:
            self.log.warn('ERROR: Failed to write %s (ccdId=%d): %s\n' % (dataId, ccdId, e))

        self.log.info("Finished writing CCD %s (ccdId=%d) on %s" % (dataId, ccdId, NODE))
