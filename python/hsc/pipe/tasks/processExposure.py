import os
import sys
import collections

import pbasf2 as pbasf
import hsc.pipe.tasks.plotSetup
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import hsc.pipe.base.matches as hscMatches
from hsc.pipe.base.argumentParser import SubaruArgumentParser
from lsst.pipe.base import Struct, CmdLineTask
from lsst.pex.config import Config, ConfigurableField
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask
from hsc.pipe.base.mpi import abortOnError, thisNode


class MpiArgumentParser(SubaruArgumentParser):
    @abortOnError
    def _makeDataRefList(self, *args, **kwargs):
        # We don't want all the MPI jobs to go reading the registry at once
        comm = pbasf.Comm()
        rank = comm.rank
        root = 0
        if rank == root:
            dataRefList = super(MpiArgumentParser, self)._makeDataRefList(*args, **kwargs)
        else:
            dataRefList = None
        if comm.size > 1:
            num = pbasf.Broadcast(len(dataRefList), root=root)
            if rank != root:
                # Ensure there's the same number of entries
                dataRefList = [None] * num
        return dataRefList

class ProcessExposureConfig(Config):
    processCcd = ConfigurableField(target=SubaruProcessCcdTask, doc="CCD processing task")

    def setDefaults(self):
        # We will do persistence ourselves
        self.processCcd.isr.doWrite = False
        self.processCcd.doWriteCalibrate = False
        self.processCcd.doWriteSources = False
        self.processCcd.doFinalWrite = False


class ProcessExposureTask(CmdLineTask):
    ConfigClass = ProcessExposureConfig
    _DefaultName = "processExposure"

    def __init__(self, **kwargs):
        super(ProcessExposureTask, self).__init__(**kwargs)
        self.comm = pbasf.Comm()
        self.rank = self.comm.rank
        self.root = 0
        self.processCcd = self.makeSubtask("processCcd")
        self.resultCache = dict() # Cache to remember results for saving

    @classmethod
    def _makeArgumentParser(cls):
        return MpiArgumentParser(name="processExposure", dataRefLevel="visit")

    @abortOnError
    def run(self, expRef):
        if self.rank == self.root:
            ccdRefList = dict([(ccdRef.get("ccdExposureId"), ccdRef) for ccdRef in expRef.subItems("ccd")
                               if ccdRef.datasetExists("raw")])
        else:
            ccdRefList = dict()

        # Scatter: process CCDs independently
        structList = pbasf.ScatterJob(comm, self.process, ccdRefList.values(), root=self.root)

        # Gathered: global WCS solution
        if self.rank == self.root:
            matchLists = self.getMatchLists(structList)
            wcsList = self.astrometricSolution(matchLists, expRef.get("camera"))
            fluxMag0 = self.photometricSolution(matchLists.values(), self.getFilterName(structList))
            ccdIdList = None
        else:
            matchLists = None
            wcsList = None
            fluxMag0 = None
            ccdIdList = self.resultsCache.keys()

        # Scatter with data from root: save CCDs with update astrometric/photometric solutions
        query = lambda ccdId: Struct(ccdId=ccdId, wcs=wcsList[ccdId], fluxMag0=calib.fluxMag0,
                                     dataRef=ccdRefList[ccdId])
        pbasf.QueryToRoot(comm, self.write, query, ccdIdList, root=self.root)


    def process(self, dataRef):
        ccdId = dataRef.get("ccdExposureId")
        print "Started processing %s (ccdId=%d) on %s" % (dataRef.dataId, ccdId, thisNode())
        try:
            results = self.processCcd.runDataRefList([dataRef])
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataRef.dataId, e))
            raise

        # Cache the results (in particular, the image)
        if len(results) != 1:
            self.resultsCache[ccdId] = None
            return None
        results = results[0]
        self.resultCache[ccdId] = results
        filterName = result.exposure.getFilter().getName()

        # Reformat the matches for MPI transfer
        matches, numMatches = None, 0
        if result.matches is not None and result.matchMeta is not None:
            matches = hscMatches.matchesToCatalog(result.matches, result.matchMeta)
            numMatches = len(matches)

        print ("Finished processing %s (ccdId=%d) on %s with %d matches" %
               (dataRef.dataId, ccdId, thisNode(), numMatches))

        return Struct(ccdId=ccdId, matches=matches, filterName=filterName)

    def getFilterName(self, structList):
        filterList = [s.filterName if s is not None else None for s in structList]
        filterSet = set(filterList)
        filterSet.discard(None) # Just in case
        if len(filterSet) != 1:
            raise RuntimeError("Multiple filters over exposure: %s" % filterSet)
        return filterSet.pop()

    def getMatchLists(self, structList):
        keyValue = [(s.ccdId, hscMatches.matchesFromCatalog(s.matches,
                                                            self.processCcd.measurement.config.slots))
                    for s in structList if s is not None]
        return collections.OrderedDict(keyValue)


    def astrometricSolution(self, matchLists, cameraGeom):
        try:
            from hsc.meas.tansip.solvetansip import SolveTansipTask
            config = SolveTansipTask.ConfigClass()
            task = SolveTansipTask(name="solvetansip", config=config)
            solvetansipIn = [task.convert(ml) if ml is not None else [] for ml in matchLists.values()]
            wcsList = task.solve(instrument, cameraGeom, solvetansipIn)
        except Exception, e:
            sys.stderr.write("WARNING: Global astrometric solution failed: %s\n" % e)
            wcsList = [None] * len(matchLists)
        return dict(zip(matchLists.keys(), wcsList))

    def photometricSolution(self, matchLists, filterName):
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

    def write(self, struct):
        ccdId = struct.ccdId
        if not ccdId in self.resultCache:
            # This node didn't process this CCD, or it failed; either way, nothing we can do
            return
        dataRef = struct.dataRef
        print "Start writing %s (ccdId=%d) on %s" % (dataRef.dataId, ccdId, thisNode())
        wcs = struct.wcs
        fluxMag0 = struct.fluxMag0

        try:
            result = self.resultCache[ccdId]
            self.processCcd.write(dataRef, result, wcs=wcs, fluxMag0=fluxMag0)
            del self.resultCache[ccdId]
        except Exception, e:
            sys.stderr.write('ERROR: Failed to write %s (ccdId=%d): %s\n' % (dataRef.dataId, ccdId, e))

        print "Finished writing CCD %s (ccdId=%d) on %s" % (dataRef.dataId, ccdId, thisNode())
