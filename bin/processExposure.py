#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal


def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)



def main(instrument, rerun, frameList):
    print "Processing inst=%s rerun=%s frames=%s" % (instrument, rerun, frameList)
    try:
        for frame in frameList:
            ProcessExposure(instrument, rerun, frame)
        return 0
    except:
        pbasf.ReportError("Total catastrophic failure processing frame %s" % (frameId))
        print "Aborting due to errors."
        mpi.COMM_WORLD.Abort(1)
        return 1


def ProcessExposure(instrument, rerun, frame):
    comm = mpi.COMM_WORLD

    if instrument == "hsc":
        numCcd = 100
    elif instrument == "suprimecam":
        numCcd = 10
    else:
        raise RuntimeError("Unknown instrument: %s" % (instrument))

    # Scatter: process CCDs independently
    worker = Worker(instrument, rerun)
    matchLists = pbasf.ScatterJob(comm, worker.process, [dict(visit=frame, ccd=ccd) for ccd in range(numCcd)],
                                  root=0)

    # Together: global WCS solution
    if comm.Get_rank() == 0:
        wcsList = pbasf.SafeCall(globalWcs, instrument, matchLists)
        if not wcsList:
            sys.stderr.write("Global astrometric solution failed!\n")
            wcsList = [None] * nCCD

    # Scatter with data from root: save CCDs with WCS
    pbasf.QueryToRoot(comm, worker.write, lambda dataId: wcsList[dataId['ccd']],
                      processor.resultCache.keys(), root=0)
    return 0
#end

class Worker(object):
    """Worker to process a CCD"""
    def __init__(self, instrument, rerun):
        if self.instrument == "hsc":
            import lsst.obs.hscSim as hscSim
            self.mapper = hscSim.HscSimMapper(rerun=rerun)
            ProcessCcdTask = hscProcessCcd.HscDc2ProcessCcdTask
            overrides = "hsc.py"
            # XXX override defaults
        elif self.instrument == "suprimecam":
            import lsst.obs.suprimcecam as suprimecam
            self.mapper = suprimecam.SuprimeCamMapper(rerun=rerun)
            ProcessCcdTask = hscProcessCcd.SuprimeCamProcessCcdTask
            overrides = "suprimecam.py"
        else:
            raise RuntimeError("Unrecognised instrument: %s" % self.instrument)
        self.butler = dafPersist.ButlerFactory(mapper=mapper).create()
        self.config = ProcessCcdTask.ConfigClass()
        self.config.load(os.path.join(os.environ['HSCPIPE_DIR'], 'config', overrides))
        self.processor = ProcessCcdTask(config=self.config)

        self.resultCache = dict() # Cache to remember results for saving
    
    def process(self, dataId):
        print "Started processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())

        # We will do persistence ourselves
        self.config.doWriteIsr = False
        self.config.doWriteCalibrate = False
        self.config.doWritePhotometry = False

        try:
            self.resultCache[dataId] = self.processor.runButler(self.butler, dataId)
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            raise

        print "Finished processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())
        return self.resultCache[dataId].matches

    def write(self, dataId, wcs):
        print "Start writing %d on %s,%d" % (dataId, os.uname()[1], os.getpid())

        try:
            result = self.resultCache[dataId]
            self.processor.write(self.butler, dataId, result, wcs)
            del self.resultCache[dataId]
        except Exception, e:
            sys.stderr.write('Failed to write %s: %s\n' % (dataId, e))

        print "Finished writing CCD %d on %s,%d" % (dataId, os.uname()[1], os.getpid())


def globalWcs(instrument, matchLists):
    import hsc.meas.tansip.doTansip as tansip
    import lsst.pex.policy as pexPolicy

    if instrument == "hsc":
        import lsst.obs.hscSim as hscSim
        mapper = hscSim.HscSimMapper()
        policyName = "hsc.paf"
    else:
        import lsst.obs.suprimecam as suprimecam
        mapper = suprimecam.SuprimecamMapper()
        policyName = "suprimecam.paf"

    policy = pexPolicy.Policy(os.path.join(os.getenv("SOLVETANSIP_DIR"), "policy", policyName))
    wcs = tansip.doTansip(matchLists, policy=policy, camera=mapper.camera)
    return tansip.getwcsList(wcs)


if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
