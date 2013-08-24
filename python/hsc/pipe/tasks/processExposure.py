import sys
import collections

import hsc.pipe.tasks.plotSetup
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCg
import hsc.pipe.base.butler as hscButler
from lsst.pipe.base import Struct
from lsst.pex.config import Config, Field, ConfigurableField
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask
from hsc.pipe.tasks.photometricSolution import PhotometricSolutionTask
from hsc.pipe.base.mpi import abortOnError, thisNode, MpiTask, MpiArgumentParser
from hsc.pipe.base.pbs import PbsCmdLineTask


class ProcessExposureConfig(Config):
    processCcd = ConfigurableField(target=SubaruProcessCcdTask, doc="CCD processing task")
    photometricSolution = ConfigurableField(target=PhotometricSolutionTask, doc="Global photometric solution")
    instrument = Field(dtype=str, default="suprimecam", doc="Instrument name, for solvetansip")
    doSolveTansip = Field(dtype=bool, default=True, doc="Run solvetansip?")

    def setDefaults(self):
        # We will do persistence ourselves
        self.processCcd.isr.doWrite = False
        self.processCcd.doWriteCalibrate = False
        self.processCcd.doWriteSources = False
        self.processCcd.doWriteHeavyFootprintsInSources = False
        self.processCcd.doFinalWrite = False


class ProcessExposureTask(PbsCmdLineTask, MpiTask):
    """Process an entire exposure at once.

    We use MPI to gather the match lists for exposure-wide astrometric and
    photometric solutions.  Note that because of this, different nodes
    see different parts of the code.
    """

    RunnerClass = hscButler.ButlerTaskRunner
    ConfigClass = ProcessExposureConfig
    _DefaultName = "processExposure"

    def __init__(self, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessExposureTask, self).__init__(**kwargs)
        self.makeSubtask("processCcd")
        self.makeSubtask("photometricSolution", schema=self.processCcd.schema)
        self.resultsCache = dict() # Cache to remember results for saving

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numCcds = sum(1 for raft in parsedCmd.butler.get("camera") for ccd in afwCg.cast_Raft(raft))
        numCycles = int(numCcds/float(numNodes*numProcs) + 0.5)
        numExps = len(cls.RunnerClass.getTargetList(parsedCmd))
        return time*numExps*numCycles

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doPbs = kwargs.pop("doPbs", False)
        parser = MpiArgumentParser(name="processExposure", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", level="visit",
                               help="data ID, e.g. --id visit=12345")
        return parser

    @abortOnError
    def run(self, expRef, butler):
        """Process a single exposure, with scatter-gather-scatter using MPI.

        All nodes execute this method, though the master and slaves have different
        routes through it.  The expRef is only a DummyDataRef on the slaves.
        """
        self.butler = butler

        if self.rank == self.root:
            dataIdList = dict([(ccdRef.get("ccdExposureId"), ccdRef.dataId)
                               for ccdRef in expRef.subItems("ccd") if ccdRef.datasetExists("raw")])
        else:
            dataIdList = dict()

        # Scatter: process CCDs independently
        import pbasf2
        structList = pbasf2.ScatterJob(self.comm, self.process, dataIdList.values(), root=self.root)
        if self.rank == self.root:
            numGood = sum(1 for s in structList if s is not None)
        else:
            numGood = 0
        if self.comm.size > 1:
            numGood = pbasf2.Broadcast(self.comm, numGood, root=self.root)
        if numGood == 0:
            return

        # Gathered: global WCS solution
        if self.rank == self.root:
            matchLists = self.getMatchLists(structList)
            wcsList = self.solveAstrometry(matchLists, self.butler.mapper.camera)
            fluxMag0 = self.solvePhotometry(matchLists.values(), self.getFilterName(structList))
            ccdIdList = self.resultsCache.keys()
        else:
            matchLists = None
            wcsList = None
            fluxMag0 = None
            ccdIdList = self.resultsCache.keys()

        # Scatter with data from root: save CCDs with update astrometric/photometric solutions
        query = lambda ccdId: Struct(wcs=wcsList[ccdId], fluxMag0=fluxMag0, dataId=dataIdList[ccdId])
        pbasf2.QueryToRoot(self.comm, self.write, query, ccdIdList, root=self.root)


    def process(self, dataId):
        """Process a single CCD and save the results for a later write.

        Only slaves execute this method.
        """
        dataRef = hscButler.getDataRef(self.butler, dataId)
        ccdId = dataRef.get("ccdExposureId")
        print "Started processing %s (ccdId=%d) on %s" % (dataId, ccdId, thisNode())
        try:
            result = self.processCcd.run(dataRef)
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            return None

        # Cache the results (in particular, the image)
        self.resultsCache[ccdId] = result
        filterName = result.exposure.getFilter().getName()

        # Reformat the matches for MPI transfer
        matches, numMatches = None, 0
        if result.matches is not None and result.matchMeta is not None:
            matches = afwTable.ReferenceMatchVector(result.matches)
            numMatches = len(matches)

        print ("Finished processing %s (ccdId=%d) on %s with %d matches" %
               (dataId, ccdId, thisNode(), numMatches))

        return Struct(ccdId=ccdId, matches=matches, filterName=filterName)

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

    def getMatchLists(self, structList):
        """Generate a list of matches for each CCD from the list of structs returned by process().

        The matches are reconsituted from the transfer format.

        Only the master executes this method, as the structList is only valid there.
        """
        keyValue = [(s.ccdId, s.matches) for s in structList if s is not None]
        return collections.OrderedDict(sorted(keyValue, key=lambda kv: kv[0]))

    def solveAstrometry(self, matchLists, cameraGeom):
        """Determine a global astrometric solution for the exposure.

        Only the master executes this method, as the matchLists is only valid there.
        """
        wcsList = [None] * len(matchLists)
        if self.config.doSolveTansip:
            try:
                from hsc.meas.tansip.solvetansip import SolveTansipTask
                config = SolveTansipTask.ConfigClass()
                task = SolveTansipTask(name="solvetansip", config=config)
                solvetansipIn = [task.convert(ml) if ml is not None else [] for ml in matchLists.values()]
                wcsList = task.solve(self.config.instrument, cameraGeom, solvetansipIn)
            except Exception, e:
                sys.stderr.write("WARNING: Global astrometric solution failed: %s\n" % e)
        else:
            self.log.info("solvetansip disabled in configuration")
        return dict(zip(matchLists.keys(), wcsList))

    def solvePhotometry(self, matchLists, filterName):
        try:
            return self.photometricSolution.run(matchLists, filterName)
        except Exception, e:
            self.log.warn("Failed to determine global photometric zero-point: %s" % e)
        return None

    def write(self, ccdId, struct):
        """Write the outputs.

        The cached results are written along with revised astrometric and photometric solutions.

        This method is only executed on the slaves.
        """
        if not ccdId in self.resultsCache:
            # This node didn't process this CCD, or it failed; either way, nothing we can do
            return
        dataId = struct.dataId
        dataRef = hscButler.getDataRef(self.butler, dataId)
        print "Start writing %s (ccdId=%d) on %s" % (dataId, ccdId, thisNode())
        wcs = struct.wcs
        fluxMag0 = struct.fluxMag0

        try:
            result = self.resultsCache[ccdId]
            self.processCcd.write(dataRef, result, wcs=wcs, fluxMag0=fluxMag0)
            del self.resultsCache[ccdId]
        except Exception, e:
            sys.stderr.write('ERROR: Failed to write %s (ccdId=%d): %s\n' % (dataId, ccdId, e))

        print "Finished writing CCD %s (ccdId=%d) on %s" % (dataId, ccdId, thisNode())
