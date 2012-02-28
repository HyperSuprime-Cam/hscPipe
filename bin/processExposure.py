#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal

import hsc.pipe.base.camera as hscCamera
import hsc.pipe.tasks.processCcd as hscProcessCcd

def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)



def main(instrument, rerun, frameList):
    print "Processing inst=%s rerun=%s frames=%s" % (instrument, rerun, frameList)
    try:
        for frame in frameList:
            print "Processing frame %d" % frame
            ProcessExposure(instrument, rerun, frame)
            print "Done processing frame %d" % frame
    except:
        pbasf.ReportError("Total catastrophic failure processing frame %s" % frame)
        print "Aborting due to errors."
        mpi.COMM_WORLD.Abort(1)
        return 1

def ProcessExposure(instrument, rerun, frame):
    comm = mpi.COMM_WORLD

    if instrument.lower() in ("hsc"):
        ProcessCcdTask = hscProcessCcd.HscDc2ProcessCcdTask
        overrides = "hsc.py"
    elif instrument == "suprimecam":
        ProcessCcdTask = hscProcessCcd.SuprimeCamProcessCcdTask
        overrides = "suprimecam.py"
    else:
        raise RuntimeError("Unknown instrument: %s" % (instrument))

    butler = hscCamera.getButler(instrument, rerun)
    dataIdList = [{'visit': frame, 'ccd': ccd} for ccd in range(hscCamera.getNumCcds(instrument))]

    config = ProcessCcdTask.ConfigClass()
    config.load(os.path.join(os.environ['HSCPIPE_DIR'], 'config', overrides))
    processor = ProcessCcdTask(config=config)

    # Scatter: process CCDs independently
    worker = Worker(butler, processor)
    matchLists = pbasf.ScatterJob(comm, worker.process, dataIdList, root=0)

    # Together: global WCS solution
    if comm.Get_rank() == 0:
        wcsList = pbasf.SafeCall(globalWcs, instrument, matchLists) if False else None
        if not wcsList:
            sys.stderr.write("WARNING: Global astrometric solution failed!\n")
            wcsList = [None] * len(dataIdList)

    # Scatter with data from root: save CCDs with WCS
    pbasf.QueryToRoot(comm, worker.write, lambda dataId: wcsList[dataId['ccd']], dataIdList, root=0)

class Worker(object):
    """Worker to process a CCD"""
    def __init__(self, butler, processor):
        self.butler = butler
        self.processor = processor
        self.resultCache = dict() # Cache to remember results for saving
    
    def process(self, dataId):
        print "Started processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())

        # ProcessCcdTask wants a dataRef, but they don't pickle so we need to reconstruct it from the dataId
        dataRefList = [ref for ref in self.butler.subset(datasetType='raw', **dataId)]
        assert len(dataRefList) == 1
        dataRef = dataRefList[0]

        # We will do persistence ourselves
        self.processor.config.doWriteIsr = False
        self.processor.config.doWriteCalibrate = False
        self.processor.config.doWritePhotometry = False
            

        try:
            self.resultCache[dataId['ccd']] = self.processor.run(dataRef)
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            raise

        print "Finished processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())
        return [m for m in self.resultCache[dataId['ccd']].matches]

    def write(self, dataId, wcs):
        if not dataId['ccd'] in self.resultCache:
            return
        print "Start writing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())

        try:
            result = self.resultCache[dataId['ccd']]
            self.processor.write(self.butler, dataId, result, wcs)
            del self.resultCache[dataId['ccd']]
        except Exception, e:
            sys.stderr.write('ERROR: Failed to write %s: %s\n' % (dataId, e))

        print "Finished writing CCD %s on %s,%d" % (dataId, os.uname()[1], os.getpid())


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
    policy.set('NCCD', len(matchLists))
    wcs = tansip.doTansip(matchLists, policy=policy, camera=mapper.camera)
    return tansip.getwcsList(wcs)


if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
