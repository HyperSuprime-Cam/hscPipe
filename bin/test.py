print "importing sys"
import sys
import mpi4py.MPI as mpi
print "importing pbasf"
import pbasf2 as pbasf
print "importing various others"
import os
import signal
import time

def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)
    
def main(instrument, rerun, lFrameId):
    print "running inst=%s rerun=%s frames=%s" % (instrument, rerun, lFrameId)
    try:
        for frameId in lFrameId:
            ProcessFrame(instrument, rerun, frameId)
        return 0
    except:
        pbasf.ReportError("Total catastrophic failure processing frame %s" % (frameId))
        print "THIS ERROR SHALL NOT HAVE APPEARED."
        mpi.COMM_WORLD.Abort(1)
        return 1

# end


def ProcessFrame(instrument, rerun, frameId):
    comm = mpi.COMM_WORLD

    if instrument == "hsc":
        nCCD = 100
    elif instrument == "suprimecam":
        nCCD = 10
    else:
        raise RuntimeError("unknown instrument: %s" % (instrument))
    
    tester = TestWorker(rerun=rerun, instrument=instrument)
    output = pbasf.ScatterJob(comm, tester, [(frameId, ccdId) for ccdId in range(nCCD)], root=0)

    print output
    return 0

class TestWorker:
    def __init__(self, rerun=None, instrument="hsc"):
        self.rerun = rerun
        self.instrument = instrument
        self.ccdIdToCache = dict()

    def __call__(self, t_frameId_ccdId):
        ccdId = t_frameId_ccdId[1]
        frameId = t_frameId_ccdId[0]
        print "Test processing %d,%d on %s,%d" % (frameId, ccdId, os.uname()[1], os.getpid())
        return (frameId, ccdId)

if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
