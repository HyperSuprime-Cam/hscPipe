#!/usr/bin/env python

# For help, invoke this script without arguments.
#
import os
import stat
import sys

def CreatePBSScript \
        (   lFrameId
            ,   sRerunName
            ,   sJobName
            ,   nNodes
            ,   nProcsPerNode
            ,   secCpuTime
            ,   sOutputDir
            ,   sQueueName = None
            ,   sInstrument = "hsc"
            ,   sProc="blammo"
            ):
    """
        Create PBS script that launches pbasf.
        Return: Path of the PBS script
    """

    sPBSScriptName   = "pbs-job-%s-%s.sh" % (sRerunName, "+".join(lFrameId))

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
    frameCpuTime = len(lFrameId) * nCcd * secCpuTime 
    frameWallTime = frameCpuTime / (float(nNodes)*nProcsPerNode)
    
    f = open(sPBSScriptName, "w")

    print >>f, "#!/bin/bash"
    print >>f, "#   Post this job with `qsub -V $0'"
    print >>f, "#PBS -l nodes=%d:ppn=%d"     % (nNodes, nProcsPerNode)
    print >>f, "#PBS -l walltime=%d"         % (2*frameWallTime,     )
    print >>f, "#PBS -o %s"                  % (sOutputDir,          )
    print >>f, "#PBS -N %s"                  % (sJobName,            )
    if sQueueName:
        print >>f, "#PBS -q %s"                  % (sQueueName           )
    print >>f, "#PBS -j oe"
    print >>f, "#PBS -W umask=02"
    print >>f, ""
    print >>f, "setup mpiexec"
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
    print >>f, "mpiexec --verbose python %s %s %s %s" % (sPBASFScriptPath,
                                                         sInstrument, sRerunName,
                                                         ' '.join(str(i) for i in lFrameId))
    print >>f, ""

    del f

    os.chmod \
    (   sPBSScriptName
    ,   stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR
    )

    return sPBSScriptName


#
# command line procs
#
g_shortopts  = ""
g_longopts   = []
g_dispacher  = dict()
g_globalVars = dict()
g_descript   = []

def CmdOption_Help(val):
    g_globalVars["help"] = True

g_shortopts +=    "h"
g_dispacher[     '-h'   ] = CmdOption_Help
g_longopts.append("help")
g_dispacher[    '--help'] = CmdOption_Help
g_globalVars[     "help"] = False
g_descript.append \
(   "\n-h, --help"
+   "\n    Show this message and exit."
)

def CmdOption_Version(val):
    g_globalVars["version"] = True

g_shortopts +=    "V"
g_dispacher[     '-V'      ] = CmdOption_Version
g_longopts.append("version")
g_dispacher[    '--version'] = CmdOption_Version
g_globalVars[     "version"] = False
g_descript.append \
(   "\n-V, --version"
+   "\n    Show version."
)


def CmdOption_RerunName(val):
    g_globalVars["rerun-name"] = val

g_shortopts +=    "r:"
g_dispacher[     '-r'          ] = CmdOption_RerunName
g_longopts.append("rerun-name=")
g_dispacher[    '--rerun-name' ] = CmdOption_RerunName
g_globalVars[     "rerun-name" ] = os.getenv('USER')
g_descript.append \
(   "\n-r, --rerun-name"
+   "\n    Specify rerun name (%s)." % g_globalVars["rerun-name"]
)


def CmdOption_JobName(val):
    g_globalVars["job-name"] = val

g_shortopts +=    "j:"
g_dispacher[     '-j'        ] = CmdOption_JobName
g_longopts.append("job-name=")
g_dispacher[    '--job-name' ] = CmdOption_JobName
g_globalVars[     "job-name" ] = "reduceFrames"
g_descript.append \
(   "\n-j, --job-name"
+   "\n    Specify PBS job name (%s)." % g_globalVars["job-name"]
)


def CmdOption_Nodes(val):
    g_globalVars["nodes"] = int(val)

g_shortopts +=    "n:"
g_dispacher[     '-n'     ] = CmdOption_Nodes
g_longopts.append("nodes=")
g_dispacher[    '--nodes' ] = CmdOption_Nodes
g_globalVars[     "nodes" ] = 4
g_descript.append \
(   "\n-n, --nodes"
+   "\n    The number of nodes on which to let the jub run (%d)." % g_globalVars["nodes"]
)


def CmdOption_nProcsPerNode(val):
    g_globalVars["procs-per-node"] = int(val)

g_shortopts +=    "p:"
g_dispacher[     '-p'              ] = CmdOption_nProcsPerNode
g_longopts.append("procs-per-node=")
g_dispacher[    '--procs-per-node' ] = CmdOption_nProcsPerNode
g_globalVars[     "procs-per-node" ] = 8
g_descript.append \
(   "\n-p, --procs-per-node"
+   "\n    Processes per node. (%d)." % g_globalVars["procs-per-node"]
)


def CmdOption_CpuTime(val):
    g_globalVars["cpu-time"] = int(val)

g_shortopts +=    "w:"
g_dispacher[     '-w'         ] = CmdOption_CpuTime
g_longopts.append("cpu-time=")
g_dispacher[    '--cpu-time' ] = CmdOption_CpuTime
g_globalVars[     "cpu-time" ] = 120
g_descript.append \
(   "\n-w, --cpu-time"
+   "\n    Roughly expected time of execution in seconds. (%d)." % g_globalVars["cpu-time"]
)


def CmdOption_OutputDir(val):
    g_globalVars["output-dir"] = val

g_shortopts +=    "o:"
g_dispacher[     '-o'          ] = CmdOption_OutputDir
g_longopts.append("output-dir=")
g_dispacher[    '--output-dir' ] = CmdOption_OutputDir
g_globalVars[     "output-dir" ] = None
g_descript.append \
(   "\n-o, --output-dir"
+   "\n    Specify whereon to write stdin/err log. (current directory)."
)

def CmdOption_Instrument(val):
    g_globalVars["instrument"] = val

g_shortopts +=    "i:"
g_dispacher[     '-i'          ] = CmdOption_Instrument
g_longopts.append("instrument=")
g_dispacher[    '--instrument' ] = CmdOption_Instrument
g_globalVars[     "instrument" ] = "hsc"
g_descript.append \
(   "\n-i, --instrument"
+   "\n    Specify which instrument to reduce for (hsc)."
)

# Help Messages
def OnHelp():
    print "Reduce frames. The job is posted to PBS queue."
    print ""
    print "Command line:"
    print "    reduceFrames [OPTIONS...] FRAMEIDS..."
    print ""
    print "FRAMEIDS... is a list of integers identifying frames."
    print "OPTIONS... are as follows."

    for i in g_descript:
        print i

def OnVersion():
    print "reduceFrames version 100828-01"


################################################
# Entry point

#
# analyze options
#

import getopt
import sys

options, args = getopt.gnu_getopt(sys.argv, g_shortopts, g_longopts)

for key, val in options:
    g_dispacher[key](val);

# if no args given, print help
if len(args) <= 1 and not g_globalVars["version"]:
    g_globalVars["help"] = True

#
# exit if help or version
#

if g_globalVars["help"] or g_globalVars["version"]:
    if g_globalVars["help"]:
        OnHelp()
        print "\n"
    OnVersion()
    sys.exit(0)

#
# argments other than options are all frame IDs.
#

lFrameId = [ str(int(i)) for i in args[1:] ]

#
# Create PBS Script
#
sScriptName = CreatePBSScript \
(   lFrameId      = lFrameId
,   sRerunName    = g_globalVars["rerun-name"]
,   sJobName      = "%s-%s" % (g_globalVars["job-name"], g_globalVars["instrument"]) 
,   nNodes        = g_globalVars["nodes"]
,   nProcsPerNode = g_globalVars["procs-per-node"]
,   secCpuTime    = g_globalVars["cpu-time"]
,   sOutputDir    = g_globalVars["output-dir"]
,   sInstrument   = g_globalVars["instrument"]
,   sProc         = "DC2-phase1.py"
)

#
# Post the script
#

os.execlp('qsub', 'qsub', '-V', sScriptName)

# eof
