#!/usr/bin/env python

# For help, invoke this script without arguments.
#
import os
import stat
import sys
import tempfile

def CreatePBSScript \
        (   sField
            ,   sFilter
            ,   sDestWcs
            ,   sRerunName
            ,   sJobName
            ,   nNodes
            ,   nProcsPerNode
            ,   secCpuTime
            ,   sOutputDir
            ,   sQueueName = None
            ,   sInstrument = "hsc"
	    ,   doMatchPsf = False
            ,   sProc="blammo"
            ):
    """
        Create PBS script that launches pbasf.
        Return: Path of the PBS script
    """

    sPBASFScriptPath = os.path.join \
    (   os.path.dirname(os.path.abspath(sys.argv[0]))
    ,   sProc
    )

    sCwd = os.getcwd()

    if not sOutputDir:
        sOutputDir = sCwd

    if sInstrument == "hsc":
        nCcd = 100
    else:
        nCcd = 10
    frameCpuTime = nCcd * secCpuTime 
    frameWallTime = frameCpuTime / (float(nNodes)*nProcsPerNode)
    
    fd, sPBSScriptName = tempfile.mkstemp()
    f = os.fdopen(fd, "w")

    print >>f, "#!/bin/bash"
    print >>f, "#   Post this job with `qsub -V $0'"
    print >>f, "#PBS -l nodes=%d:ppn=%d"     % (nNodes, nProcsPerNode)
    print >>f, "#PBS -l walltime=%d"         % (3*frameWallTime,     )
    print >>f, "#PBS -o %s"                  % (sOutputDir,          )
    print >>f, "#PBS -N %s"                  % (sJobName,            )
    if sQueueName:
        print >>f, "#PBS -q %s"                  % (sQueueName           )
    print >>f, "#PBS -j oe"
    print >>f, "#PBS -W umask=02"
    print >>f, ""
    print >>f, "echo \"mpiexec is at: $(which mpiexec)\""
    print >>f, ""
    print >>f, "ulimit -a"
    print >>f, "umask 02"
    print >>f, "echo 'umask: ' $(umask)"
    print >>f, ""
    print >>f, "eups list -s"
    print >>f, ""
    print >>f, "export"
    print >>f, ""
    print >>f, "cd %s" % (sCwd,)
    print >>f, ""

    matchPsf = "--doMatchPsf" if doMatchPsf else ""
    if destWcs:
        print >>f, "mpiexec --verbose python %s --instrument=%s --rerun=%s --program=%s --filter=%s --workDirRoot=%s %s --destWcs=%s" % (sPBASFScriptPath,
                                                sInstrument, sRerunName,
                                                sField, sFilter, sOutputDir, matchPsf, sDestWcs)
    else:
        print >>f, "mpiexec --verbose python %s --instrument=%s --rerun=%s --program=%s --filter=%s --workDirRoot=%s %s" % (sPBASFScriptPath,
                                   sInstrument, sRerunName,
                                   sField, sFilter, sOutputDir, matchPsf)
    print >>f, ""

    del f

    os.chmod \
    (   sPBSScriptName
    ,   stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
    )

    return sPBSScriptName


################################################
# Entry point

#
# analyze options
#

import optparse
import sys

usage = "usage: %prog [options] field filter [destWcs]"
parser = optparse.OptionParser(usage=usage)
parser.add_option("-q", "--queue",
                  type=str, default=None,
                  help="Specify PBS queue name (%default).")
parser.add_option("-j", "--job",
                  type=str, default="reduceStack",
                  help="Specify PBS job name (%default).")
parser.add_option("-n", "--nodes",
                  type=int, default=4,
                  help="The number of nodes on which to let the jub run (%default).")
parser.add_option("-p", "--procs-per-node",
                  type=int, default=8,
                  help="Processes per node (%default).")
parser.add_option("-w", "--cpu-time",
                  type=int, default=120,
                  help="Roughly expected time of execution in seconds (%default).")
parser.add_option("-o", "--output-dir",
                  type=str, default=".",
                  help="Specify whereon to write stdin/err log (%default).")
parser.add_option("-r", "--rerun",
                  type=str, default=os.getenv('USER'),
                  help="Specify rerun name (%default).")
parser.add_option("-i", "--instrument",
                  type=str, default="hsc",
                  help="Specify which instrument to reduce for (%default).")
parser.add_option("-m", "--doMatchPsf",
		  default=False, action='store_true',
		  help="match PSFs before stacking (default=%default)")

(opts, args) = parser.parse_args()

if (len(args) < 2):
    parser.print_help()
    sys.exit(1)

field = args[0]
filter = args[1]
if (len(args) == 3):
    destWcs = args[2]
else:
    destWcs = None

#
# Create PBS Script
#
sScriptName = CreatePBSScript \
(   sField        = field
,   sFilter       = filter
,   sDestWcs      = destWcs
,   sRerunName    = opts.rerun
,   sJobName      = "%s-%s" % (opts.job, opts.instrument) 
,   nNodes        = opts.nodes
,   nProcsPerNode = opts.procs_per_node
,   secCpuTime    = opts.cpu_time
,   sOutputDir    = opts.output_dir
,   sQueueName    = opts.queue
,   sInstrument   = opts.instrument
,   doMatchPsf    = opts.doMatchPsf
,   sProc         = "stackExposures.py"
)

#
# Post the script
#

os.execlp('qsub', 'qsub', '-V', sScriptName)

# eof
