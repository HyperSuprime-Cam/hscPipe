import sys
import math
import collections

import hsc.pipe.tasks.plotSetup
import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCg
import hsc.pipe.base.butler as hscButler
from lsst.pipe.base import Struct, ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask
from hsc.pipe.tasks.photometricSolution import PhotometricSolutionTask
from hsc.pipe.tasks.focusTask import ProcessFocusTask
from hsc.pipe.base.pool import abortOnError, NODE, Pool, Debugger
from hsc.pipe.base.pbs import PbsPoolTask
from hsc.meas.tansip.solvetansip import SolveTansipTask

Debugger().enabled = True


class ProcessExposureConfig(Config):
    processCcd = ConfigurableField(target=SubaruProcessCcdTask, doc="CCD processing task")
    focus = ConfigurableField(target=ProcessFocusTask, doc="Focus processing task")
    photometricSolution = ConfigurableField(target=PhotometricSolutionTask, doc="Global photometric solution")
    instrument = Field(dtype=str, default="suprimecam", doc="Instrument name, for solvetansip")
    doSolveTansip = Field(dtype=bool, default=True, doc="Run solvetansip?")
    solveTansip = ConfigurableField(target=SolveTansipTask, doc="Global astrometric solution")
    doPhotometricSolution = Field(dtype=bool, default=False, doc="Run global photometric solution?")
    doFocus = Field(dtype=bool, default=False, doc="Run focus analysis?")

    def setDefaults(self):
        # We will do persistence ourselves
        self.processCcd.isr.doWrite = False
        self.processCcd.doWriteCalibrate = False
        self.processCcd.doWriteSources = False
        self.processCcd.doWriteHeavyFootprintsInSources = False
        self.processCcd.doFinalWrite = False


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
        self.makeSubtask("focus")
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

        # Gathered: global solutions (photometry, astrometry, focus)
        matchLists = self.getMatchLists(dataIdList, structList)
        wcsList = self.solveAstrometry(matchLists, butler.mapper.camera)
        filterName = self.getFilterName(structList)
        fluxMag0 = self.solvePhotometry(matchLists.values(), filterName)
        focusMd = self.solveFocus(expRef, [s.focus for s in structList if
                                           s is not None and s.focus is not None])
        ccdIdList = dataIdList.keys()

        # Scatter with data from root: save CCDs with update astrometric/photometric solutions
        solutionList = [Struct(ccdId=ccdId, fluxMag0=fluxMag0, dataId=dataIdList[ccdId],
                               wcs=wcsList[ccdId] if ccdId in wcsList else None)
                        for ccdId in ccdIdList]
        pool.mapToPrevious(self.write, solutionList, focusMd)

    def process(self, cache, dataId):
        """Process a single CCD and save the results for a later write.

        Only slaves execute this method.
        """
        dataRef = hscButler.getDataRef(cache.butler, dataId)
        ccdId = dataRef.get("ccdExposureId")
        with self.logOperation("processing %s (ccdId=%d)" % (dataId, ccdId)):
            matches = None
            filterName = None
            focus = None
            result = None
            cache.result = None
            try:
                if self.focus.isFocus(dataRef):
                    if self.config.doFocus:
                        focus = self.focus.run(dataRef)
                else:
                    result = self.processCcd.run(dataRef)
            except Exception, e:
                self.log.warn("Failed to process %s: %s\n" % (dataId, e))
                import traceback
                traceback.print_exc()
                return None

            if result is not None:
                # Cache the results (in particular, the image)
                cache.result = result
                filterName = result.exposure.getFilter().getName()

                # Reformat the matches for MPI transfer
                if result.matches is not None and result.matchMeta is not None:
                    matches = afwTable.ReferenceMatchVector(result.matches)

            return Struct(ccdId=ccdId, matches=matches, filterName=filterName, focus=focus)

    def getFilterName(self, structList):
        """Determine the filter name from the list of structs returned by process().

        Only the master executes this method, as the structList is only valid there.
        """
        filterList = [s.filterName if s is not None else None for s in structList]
        filterSet = set(filterList)
        filterSet.discard(None) # Just in case
        if len(filterSet) == 0:
            self.log.warn("Unable to determine filter name")
            return None
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

    def solveFocus(self, dataRef, focusList):
        if not self.config.doFocus:
            return None

        camera = dataRef.get("camera")
        plotName = dataRef.get("focusPlot_filename")
        focus = self.focus.measureFocus(focusList, camera, plotName)
        self.log.info("Focus solution: %s" % (focus,))

        metadata = dafBase.PropertyList()
        metadata.set("FOCUS_CORRECTION_VALUE", focus[0])
        metadata.set("FOCUS_CORRECTION_ERROR", focus[1])
        metadata.set("FOCUS_CORRECTION_UNCORRECTED_VALUE", focus[2])
        metadata.set("FOCUS_CORRECTION_UNCORRECTED_ERROR", focus[3])

        return metadata


    def write(self, cache, struct, focusMd=None):
        """Write the outputs.

        The cached results are written along with revised astrometric and photometric solutions.

        This method is only executed on the slaves.
        """
        if not cache.result or not struct:
            # Processing must have failed: nothing we can do
            return

        dataId = struct.dataId
        ccdId = struct.ccdId
        with self.logOperation("writing CCD %s (ccdId=%d)" % (dataId, ccdId)):
            dataRef = hscButler.getDataRef(cache.butler, dataId)
            wcs = struct.wcs
            fluxMag0 = struct.fluxMag0

            self.processCcd.write(dataRef, cache.result, wcs=wcs, fluxMag0=fluxMag0, focusMd=focusMd)
            del cache.result
