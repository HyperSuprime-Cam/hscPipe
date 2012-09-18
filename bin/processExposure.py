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
import hsc.pipe.base.camera as hscCamera
import hsc.pipe.base.butler as hscButler
from hsc.pipe.tasks.processCcd import SubaruProcessCcdTask

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
    comm = mpi.COMM_WORLD

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
    matchLists = pbasf.ScatterJob(comm, worker.process, dataIdList, root=0)

    # Together: global WCS solution
    if comm.Get_rank() == 0:
        wcsList = pbasf.SafeCall(globalWcs, instrument, butler.mapper.camera, matchLists)
        if not wcsList:
            sys.stderr.write("WARNING: Global astrometric solution failed!\n")
            wcsList = [None] * len(dataIdList)

    # Scatter with data from root: save CCDs with WCS
    pbasf.QueryToRoot(comm, worker.write, lambda dataId: wcsList[dataId['ccd']], dataIdList, root=0)


Match = collections.namedtuple("Match", ["id", "ra", "dec", "x", "y", "xErr", "yErr", "flux"])


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

        try:
            dataRef = hscButler.getDataRef(self.butler, dataId)
            results = self.processor.runDataRefList([dataRef])
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            raise

        if len(results) != 1:
            self.resultCache[ccd] = None
            return []

        result = results[0]
        self.resultCache[ccd] = result

        try:
            matches = [Match(m.second.getId(), m.first.getRa().asDegrees(), m.first.getDec().asDegrees(),
                             m.second.getX(), m.second.getY(),
                             m.second.get(m.second.getTable().getCentroidErrKey()[0,0]),
                             m.second.get(m.second.getTable().getCentroidErrKey()[1,1]),
                             m.second.getPsfFlux()) for m in result.calib.matches]
        except:
            matches = []

        print "Finished processing %s on %s with %d matches" % (dataId, thisNode(), len(matches))

        return matches

    def write(self, dataId, wcs):
        if not dataId['ccd'] in self.resultCache:
            # This node didn't process this CCD, or it failed; either way, nothing we can do
            return
        print "Start writing %s on %s" % (dataId, thisNode())

        try:
            result = self.resultCache[dataId['ccd']]
            dataRef = hscButler.getDataRef(self.butler, dataId)
            self.processor.write(dataRef, result, wcs)
            del self.resultCache[dataId['ccd']]
        except Exception, e:
            sys.stderr.write('ERROR: Failed to write %s: %s\n' % (dataId, e))

        print "Finished writing CCD %s on %s" % (dataId, thisNode())


def globalWcs(instrument, cameraGeom, matchLists):
    import hsc.meas.tansip as tansip
    from hsc.meas.tansip.solvetansip import SolveTansipTask
    config = SolveTansipTask.ConfigClass()
    task = SolveTansipTask(name="solvetansip", config=config)

    def solveTansipTranslate(matchList):
        if matchList is None:
            return []
        return [tansip.SourceMatch(m.id, afwCoord.IcrsCoord(afwGeom.Angle(m.ra, afwGeom.degrees),
                                                            afwGeom.Angle(m.dec, afwGeom.degrees)),
                                   afwGeom.Point2D(m.x, m.y), afwGeom.Point2D(m.xErr, m.yErr), m.flux)
                for m in matchList]

    solvetansipIn = [solveTansipTranslate(ml) for ml in matchLists]

    return task.solve(instrument, cameraGeom, solvetansipIn)

if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
