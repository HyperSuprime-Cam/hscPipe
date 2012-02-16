#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal

import hsc.pipe.tasks.processCcd as hscProcessCcd
import lsst.daf.persistence as dafPersist

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
        pbasf.ReportError("Total catastrophic failure processing frame %s" % frame)
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

    if instrument == "hsc":
        import lsst.obs.hscSim as hscSim
        mapper = hscSim.HscSimMapper(rerun=rerun)
        ProcessCcdTask = hscProcessCcd.HscDc2ProcessCcdTask
        overrides = "hsc.py"
    elif instrument == "suprimecam":
        import lsst.obs.suprimecam as suprimecam
        mapper = suprimecam.SuprimecamMapper(rerun=rerun)
        ProcessCcdTask = hscProcessCcd.SuprimeCamProcessCcdTask
        overrides = "suprimecam.py"
    else:
        raise RuntimeError("Unknown instrument: %s" % (instrument))

    butler = dafPersist.ButlerFactory(mapper=mapper).create()
    config = ProcessCcdTask.ConfigClass()
    config.load(os.path.join(os.environ['HSCPIPE_DIR'], 'config', overrides))
    processor = ProcessCcdTask(config=config)

    # Scatter: process CCDs independently
    worker = Worker(processor)
    dataRefs = [ref for ref in butler.subset(datasetType='raw', visit=frame)]
    dataRefs.sort(key=lambda ref: ref.dataId['ccd']) # Ensure data references are in CCD order

    # XXX Exclude rotated CCDs for now
    if instrument == "hsc":
        dataRefs = [ref for ref in dataRefs if ref.dataId['ccd'] < 100]
    
    matchLists = pbasf.ScatterJob(comm, worker.process, dataRefs, root=0)

    # Together: global WCS solution
    if comm.Get_rank() == 0:
        wcsList = pbasf.SafeCall(globalWcs, instrument, matchLists)
        if not wcsList:
            sys.stderr.write("Global astrometric solution failed!\n")
            wcsList = [None] * len(dataRefs)

    # Scatter with data from root: save CCDs with WCS
    pbasf.QueryToRoot(comm, worker.write, lambda dataRef: wcsList[dataRef.dataId['ccd']],
                      worker.resultCache.keys(), root=0)
    return 0
#end

class Worker(object):
    """Worker to process a CCD"""
    def __init__(self, processor):
        self.processor = processor
        self.resultCache = dict() # Cache to remember results for saving
    
    def process(self, dataRef):
        dataId = dataRef.dataid
        print "Started processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())

        # We will do persistence ourselves
        self.processor.config.doWriteIsr = False
        self.processor.config.doWriteCalibrate = False
        self.processor.config.doWritePhotometry = False

        try:
            self.resultCache[dataId] = self.processor.run(dataRef)
        except Exception, e:
            sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
            raise

        print "Finished processing %s on %s,%d" % (dataId, os.uname()[1], os.getpid())
        return self.resultCache[dataId].matches

    def write(self, dataRef, wcs):
        dataId = dataRef.dataId
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
