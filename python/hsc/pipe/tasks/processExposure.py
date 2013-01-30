import os
import sys
import collections

import pbasf2 as pbasf
import hsc.pipe.tasks.plotSetup
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import hsc.pipe.base.matches as hscMatches
import hsc.pipe.base.butler as hscButler
from lsst.pipe.base import Struct, CmdLineTask
from lsst.pex.config import Config, Field, ConfigurableField
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask
from hsc.pipe.base.mpi import abortOnError, thisNode, MpiTask, MpiArgumentParser


class ProcessExposureConfig(Config):
    processCcd = ConfigurableField(target=SubaruProcessCcdTask, doc="CCD processing task")
    instrument = Field(dtype=str, default="suprimecam", doc="Instrument name, for solvetansip")

    def setDefaults(self):
        # We will do persistence ourselves
        self.processCcd.isr.doWrite = False
        self.processCcd.doWriteCalibrate = False
        self.processCcd.doWriteSources = False
        self.processCcd.doWriteHeavyFootprintsInSources = False
        self.processCcd.doFinalWrite = False


class ProcessExposureTask(MpiTask):
    """Process an entire exposure at once.

    We use MPI to gather the match lists for exposure-wide astrometric and
    photometric solutions.  Note that because of this, different nodes
    see different parts of the code.
    """

    ConfigClass = ProcessExposureConfig
    _DefaultName = "processExposure"

    def __init__(self, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessExposureTask, self).__init__(**kwargs)
        self.makeSubtask("processCcd")
        self.resultsCache = dict() # Cache to remember results for saving

    def runDataRefList(self, *args, **kwargs):
        """Save the butler.

        All nodes execute this method.
        """
        self.butler = self.parsedCmd.butler
        super(ProcessExposureTask, self).runDataRefList(*args, **kwargs)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        return MpiArgumentParser(name="processExposure", dataRefLevel="visit", *args, **kwargs)

    @abortOnError
    def run(self, expRef):
        """Process a single exposure, with scatter-gather-scatter using MPI.

        All nodes execute this method, though the master and slaves have different
        routes through it.  The expRef is only a DummyDataRef on the slaves.
        """

        if self.rank == self.root:
            dataIdList = dict([(ccdRef.get("ccdExposureId"), ccdRef.dataId)
                               for ccdRef in expRef.subItems("ccd") if ccdRef.datasetExists("raw")])
        else:
            dataIdList = dict()

        # Scatter: process CCDs independently
        structList = pbasf.ScatterJob(self.comm, self.process, dataIdList.values(), root=self.root)
        if self.rank == self.root:
            numGood = sum(1 for s in structList if s is not None)
        else:
            numGood = 0
        if self.comm.size > 1:
            numGood = pbasf.Broadcast(self.comm, numGood, root=self.root)
        if numGood == 0:
            return

        # Gathered: global WCS solution
        if self.rank == self.root:
            matchLists = self.getMatchLists(structList)
            wcsList = self.astrometricSolution(matchLists, self.butler.mapper.camera)
            fluxMag0 = self.photometricSolution(matchLists.values(), self.getFilterName(structList))
            ccdIdList = None
        else:
            matchLists = None
            wcsList = None
            fluxMag0 = None
            ccdIdList = self.resultsCache.keys()

        # Scatter with data from root: save CCDs with update astrometric/photometric solutions
        query = lambda ccdId: Struct(wcs=wcsList[ccdId], fluxMag0=fluxMag0, dataId=dataIdList[ccdId])
        pbasf.QueryToRoot(self.comm, self.write, query, ccdIdList, root=self.root)


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
            raise

        # Cache the results (in particular, the image)
        self.resultsCache[ccdId] = result
        filterName = result.exposure.getFilter().getName()

        # Reformat the matches for MPI transfer
        matches, numMatches = None, 0
        if result.matches is not None and result.matchMeta is not None:
            matches = hscMatches.matchesToCatalog(result.matches, result.matchMeta)
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
        keyValue = [(s.ccdId, hscMatches.matchesFromCatalog(s.matches,
                                                            self.processCcd.measurement.config.slots))
                    for s in structList if s is not None]
        return collections.OrderedDict(sorted(keyValue, key=lambda kv: kv[0]))

    def astrometricSolution(self, matchLists, cameraGeom):
        """Determine a global astrometric solution for the exposure.

        Only the master executes this method, as the matchLists is only valid there.
        """
        try:
            from hsc.meas.tansip.solvetansip import SolveTansipTask
            config = SolveTansipTask.ConfigClass()
            task = SolveTansipTask(name="solvetansip", config=config)
            solvetansipIn = [task.convert(ml) if ml is not None else [] for ml in matchLists.values()]
            wcsList = task.solve(self.config.instrument, cameraGeom, solvetansipIn)
        except Exception, e:
            sys.stderr.write("WARNING: Global astrometric solution failed: %s\n" % e)
            wcsList = [None] * len(matchLists)
        return dict(zip(matchLists.keys(), wcsList))

    def photometricSolution(self, matchLists, filterName):
        """Determine a global photometric solution for the exposure.

        The current implementation simply runs the general 'photocal' to get a single zero-point.

        Only the master executes this method, as the matchLists is only valid there.
        """
        photocal = self.processCcd.calibrate.photocal

        # Ensure all the match lists have the same schema, by copying them into the same table
        template = None
        for ml in matchLists:
            if ml is not None and len(ml) > 0:
                template = ml[0]
                break
        if template is None:
            raise RuntimeError("No matches provided")
        ref = afwTable.SimpleTable.make(template.first.schema)
        src = afwTable.SourceTable.make(template.second.schema)
        for prop in ("Centroid", "Shape", "PsfFlux", "ApFlux", "ModelFlux", "InstFlux"):
            getter = getattr(template.second.table, "get" + prop + "Definition")
            setter = getattr(src, "define" + prop)
            setter(getter())

        refNames = ref.schema.getNames()
        srcNames = src.schema.getNames()
        refNewKeys = [ref.schema.find(k).key for k in refNames]
        srcNewKeys = [src.schema.find(k).key for k in srcNames]
        matches = []
        for i, ml in enumerate(matchLists):
            if ml is None or len(ml) == 0:
                continue
            try:
                refOldKeys = [ml[0].first.schema.find(k).key for k in refNames]
                srcOldKeys = [ml[0].second.schema.find(k).key for k in srcNames]
            except:
                # Something's wrong with the schema; punt
                sys.stderr.write("Error with schema on matchlist %d: ignoring %d matches\n" % (i, len(ml)))
                continue

            for m in ml:
                newRef = ref.makeRecord()
                for old, new in zip(refOldKeys, refNewKeys):
                    newRef.set(new, m.first.get(old))
                newSrc = src.makeRecord()
                for old, new in zip(srcOldKeys, srcNewKeys):
                    newSrc.set(new, m.second.get(old))
                matches.append(afwTable.ReferenceMatch(newRef, newSrc, m.distance))

        try:
            class DummyExposure(object):
                """Quacks like an lsst.afw.image.Exposure, for the purposes of PhotoCal."""
                def __init__(self, filterName):
                    self._filterName = filterName
                def getFilter(self):
                    return afwImage.Filter(filterName)

            result = photocal.run(DummyExposure(filterName), matches)
        except Exception, e:
            self.log.warn("Failed to determine global photometric zero-point: %s" % e)
            return None

        self.log.info("Global photometric zero-point: %f" % result.calib.getMagnitude(1.0))
        return result.calib.getFluxMag0()

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
