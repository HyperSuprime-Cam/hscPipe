import sys
import boostmpi
import MPIFlowCtrl as flow
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
        flow.ReportError("Total catastrophic failure processing frame %s" % (frameId))
        print "THIS ERROR SHALL NOT HAVE APPEARED."
        boostmpi.abort(1)
        return 1

# end


def ProcessFrame(instrument, rerun, frameId):
    comm = boostmpi.world

    if instrument == "hsc":
        nCCD = 100
    elif instrument == "suprimecam":
        nCCD = 10
    else:
        raise RuntimeError("unknown instrument: %s" % (instrument))
    
    # create ccdId's
    lCcdId = range(nCCD)

    # phase 1
    phase1 = Phase1Worker(rerun=rerun, instrument=instrument)
    matchListAll = flow.ScatterJob \
    (   comm
    ,   phase1                                 # worker
    ,   [(frameId, ccdId) for ccdId in lCcdId] # input list
    ,   root=0                                 # domina
    )

    # phase 2
    if boostmpi.rank == 0:
        resultWcs = flow.SafeCall(phase2, instrument, matchListAll)
        if not resultWcs:
            sys.stderr.write("no global astrometric solution!!\n")
            resultWcs = [None] * nCCD
        ccdIdToWcs = dict(zip(lCcdId, resultWcs))

    # phase 3
    flow.QueryToRoot \
    (   comm
    ,   Phase3Worker(phase1.ccdIdToCache) # worker
    ,   lambda ccdId: ccdIdToWcs[ccdId]   # database
    ,   phase1.ccdIdToCache.keys()        # input list
    ,   root=0                            # who'll be database
    )
    return 0
#end

class Phase1Worker:
    def __init__(self, rerun=None, instrument="hsc"):
        self.rerun = rerun
        self.instrument = instrument
        self.ccdIdToCache = dict()

    def __call__(self, t_frameId_ccdId):
        import lsst.pipette.runHsc

        ccdId = t_frameId_ccdId[1]
        frameId = t_frameId_ccdId[0]
        # process ccdId

        print "Processing %d,%d on %s" % (frameId, ccdId, os.uname()[1])

        obj = lsst.pipette.runHsc.doRun \
              (   rerun          = self.rerun
              ,   instrument     = self.instrument
              ,   frameId        = frameId
              ,   ccdId          = ccdId
              ,   doMerge        = False
              )

        # remember the obj for phase3
        self.ccdIdToCache[ccdId] = obj

        # return matchList to the root
        return obj.matchlist

# end def


def phase2(instrument, matchListAllCcd):
    from  hsc.meas.tansip.doTansip import doTansip
    import lsst.pipette.config as pipConfig

    if instrument == "hsc":
        import lsst.obs.hscSim as hscSim
        mapper = hscSim.HscSimMapper()
        policyName = "hsc.paf"
    else:
        import lsst.obs.suprimecam as suprimecam
        mapper = suprimecam.SuprimecamMapper()
        policyName = "suprimecam.paf"

    policyPath = os.path.join(os.getenv("PIPETTE_DIR"), "policy", policyName)
    fullPolicy = pipConfig.configuration(policyPath)
    policy = fullPolicy['instrumentExtras']['solveTansip'].getPolicy()
    
    #sys.stderr.write("phase 2 being sent %d matches, #0 is %s\n" % (len(matchListAllCcd), matchListAllCcd[0]))
    return doTansip \
           (   matchListAllCcd
               ,   policy=policy
               ,   camera=mapper.camera
               )
#end def


class Phase3Worker:
    def __init__(self, ccdIdToCache):
        self.ccdIdToCache = ccdIdToCache

    def __call__(self, ccdId, wcs):
        import lsst.pipette.runHsc as runHsc

        # process it
        try:
            print "Writing CCD %d on %s" % (ccdId, os.uname()[1])
            runHsc.doMergeWcs(self.ccdIdToCache[ccdId], wcs)
            del self.ccdIdToCache[ccdId]
        except Exception, e:
            sys.stderr.write('phase3 failed to merge for %s: %s\n' % (ccdId, e))
            
#end

if __name__ == "__main__":
    print "argv=", sys.argv
    instrument = sys.argv[1]
    rerun = sys.argv[2]
    frames = [int(f) for f in sys.argv[3:]]
    main(instrument, rerun, frames)
