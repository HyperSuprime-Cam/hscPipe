import sys
import boostmpi
import MPIFlowCtrl as flow
import os
import signal
import time
import subprocess

def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)
    
def main(baseDir):
    comm = boostmpi.world
    os.chdir(baseDir)

    scripts = glob.glob('qqq[0-9]*sh')
    print "running in %s: %s" % (baseDir, scripts)
    raise SystemExit()

    runner = ScriptWorker()
    try:
        ret = flow.ScatterJob \
              (   comm
                  ,   runner                                 # worker
                  ,   scripts
                  ,   root=0                                 # domina
                  )
    except:
        flow.error("Total catastrophic failure")
        print "THIS ERROR SHALL NOT HAVE APPEARED."
        boostmpi.abort(1)
        return 1

class ScriptWorker(object):
        def __call__(self, scriptName):
            ret = subprocess.check_call([scriptName])
            return ret

if __name__ == "__main__":
    print "argv=", sys.argv
    baseDir = sys.argv[1]

    main(baseDir)
