#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal
import collections

import hsc.pipe.tasks.plotSetup
import lsst.afw.geom as afwGeom
import lsst.afw.coord as afwCoord
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import hsc.pipe.base.camera as hscCamera
import hsc.pipe.base.butler as hscButler
import hsc.pipe.base.matches as hscMatches
from lsst.pipe.base import Struct
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask

os.umask(002)

def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)

def thisNode():
    return "%s,%d" % (os.uname()[1], os.getpid())


def main(instrument, rerun, frameList):
    print "Processing inst=%s rerun=%s node=%s frames=%s" % (instrument, rerun, thisNode(), frameList)
    try:
        for frame in frameList:
            print "Processing frame %d on %s" % (frame, thisNode())
            ProcessExposure(instrument, rerun, frame)
            print "Done processing frame %d on %s" % (frame, thisNode())
    except:
        pbasf.ReportError("Total catastrophic failure processing frame %s on %s" % (frame, thisNode()))
        print "Aborting due to errors."
        mpi.COMM_WORLD.Abort(1)
        return 1
    print "Done on %s" % thisNode()

def ProcessExposure(instrument, rerun, frame):
    comm = pbasf.Comm()

    # We don't need camera-specific ProcessCcdTasks anymore.  But we may want to use
    # one of the onsite variants; in that case we should try to pass an extra argument
    # and switching which task we're using here, rather than by making a new processExposure
    # script specifically for that purpose.
    ProcessCcdTask = SubaruProcessCcdTask

    butler = hscCamera.getButler(instrument, rerun=rerun)
    dataIdList = [{'visit': frame, 'ccd': ccd} for ccd in range(hscCamera.getNumCcds(instrument))]

    # FIXME: should really rely on pipe_base to do this.
    config = ProcessCcdTask.ConfigClass()
    config.load(os.path.join(os.environ['OBS_SUBARU_DIR'], 'config',
                             ProcessCcdTask._DefaultName + ".py"))
    config.load(os.path.join(os.environ['OBS_SUBARU_DIR'], 'config', 
                             instrument, ProcessCcdTask._DefaultName + ".py"))
    processor = ProcessCcdTask(config=config)

    # Scatter: process CCDs independently
    worker = Worker(butler, processor)
    structList = pbasf.ScatterJob(comm, worker.process, dataIdList, root=0)

    # Together: global WCS solution
    if comm.Get_rank() == 0:
        matchLists = [hscMatches.matchesFromCatalog(s.matches, processor.measurement.config.slots)
                      if s is not None else None for s in structList]
        filterList = [s.filterName if s is not None else None for s in structList]

        wcsList = pbasf.SafeCall(globalWcs, instrument, butler.mapper.camera, matchLists)
        if not wcsList:
            sys.stderr.write("WARNING: Global astrometric solution failed!\n")
            wcsList = [None] * len(dataIdList)
        fluxMag0 = pbasf.SafeCall(globalZeroPoint, processor, filterList, matchLists)

    # Scatter with data from root: save CCDs with WCS
    query = lambda dataId: Struct(wcs=wcsList[dataId['ccd']], fluxMag0=fluxMag0)
    pbasf.QueryToRoot(comm, worker.write, query, dataIdList, root=0)

Match = collections.namedtuple("Match", ["id", "ra", "dec", "x", "y", "xErr", "yErr", "xyFlag", "measFlux",
                                         "measFluxErr", "measFluxFlag", "refFlux"])


class Worker(object):
    """Worker to process a CCD"""
    def __init__(self, butler, processor):
        self.butler = butler
        self.processor = processor
        self.resultCache = dict() # Cache to remember results for saving
    
    def process(self, dataId):
        print "Started processing %s on %s" % (dataId, thisNode())
        ccd = dataId['ccd']

        # We will do persistence ourselves
        self.processor.config.isr.doWrite = False
        self.processor.config.doWriteCalibrate = False
        self.processor.config.doWriteSources = False
        self.processor.config.doFinalWrite = False

        try:
            dataRef = hscButler.getDataRef(self.butler, dataId)
            results = self.processor.runDataRefList([dataRef])
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            raise

        if len(results) != 1:
            self.resultCache[ccd] = None
            return None

        result = results[0]
        self.resultCache[ccd] = result
        filterName = result.exposure.getFilter().getName()

        matches, numMatches = None, 0
        if result.matches is not None and result.matchMeta is not None:
            matches = hscMatches.matchesToCatalog(result.matches, result.matchMeta)
            numMatches = len(matches)

        print "Finished processing %s on %s with %d matches" % (dataId, thisNode(), numMatches)

        return Struct(matches=matches, filterName=filterName)

    def write(self, dataId, struct):
        if not dataId['ccd'] in self.resultCache:
            # This node didn't process this CCD, or it failed; either way, nothing we can do
            return
        print "Start writing %s on %s" % (dataId, thisNode())
        wcs = struct.wcs
        fluxMag0 = struct.fluxMag0

        try:
            result = self.resultCache[dataId['ccd']]
            dataRef = hscButler.getDataRef(self.butler, dataId)
            self.processor.write(dataRef, result, wcs=wcs, fluxMag0=fluxMag0)
            del self.resultCache[dataId['ccd']]
        except Exception, e:
            sys.stderr.write('ERROR: Failed to write %s: %s\n' % (dataId, e))

        print "Finished writing CCD %s on %s" % (dataId, thisNode())


def globalWcs(instrument, cameraGeom, matchLists):
    import hsc.meas.tansip as tansip
    from hsc.meas.tansip.solvetansip import SolveTansipTask
    config = SolveTansipTask.ConfigClass()
    task = SolveTansipTask(name="solvetansip", config=config)

    solvetansipIn = [task.convert(ml) if ml is not None else [] for ml in matchLists]
    return task.solve(instrument, cameraGeom, solvetansipIn)

def globalZeroPoint(processor, filterList, matchLists):
    photocal = processor.calibrate.photocal

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
        filterSet = set(filterList)
        filterSet.discard(None) # Just in case
        if len(filterSet) != 1:
            raise RuntimeError("Multiple filters over exposure: %s" % filterSet)
        filterName = filterSet.pop()

        class DummyExposure(object):
            """Quacks like an lsst.afw.image.Exposure, for the purposes of PhotoCal."""
            def __init__(self, filterName):
                self._filterName = filterName
            def getFilter(self):
                return afwImage.Filter(filterName)

        result = photocal.run(DummyExposure(filterName), matches)
    except Exception, e:
        processor.log.warn("Failed to determine global photometric zero-point: %s" % e)
        return None

    processor.log.info("Global photometric zero-point: %f" % result.calib.getMagnitude(1.0))
    return result.calib.getFluxMag0()



if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
