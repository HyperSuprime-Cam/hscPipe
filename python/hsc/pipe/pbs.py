#!/usr/bin/env python

import os
import os.path
import stat
import sys
import tempfile

class Pbs(object):
    def __init__(self, outputDir, numNodes=1, numProcsPerNode=1, queue=None, jobName=None, wallTime=None):
        self.outputDir = outputDir
        self.numNodes = numNodes
        self.numProcsPerNode = numProcsPerNode
        self.queue = queue
        self.jobName = jobName
        self.wallTime = wallTime

    def create(command, wallTime=None, numNodes=None, numProcsPerNode=None, jobName=None):
        if wallTime is None:
            wallTime = self.wallTime
        if numNodes is None:
            numNodes = self.numNodes
        if numProcsPerNode is None:
            numProcsPerNode = self.numProcsPerNode
        if jobName is None:
            jobName = self.jobName

        fd, script = tempfile.mkstemp()
        f = os.fdopen(fd, "w")

        assert numNodes is not None and numProcsPerNode is not None
        if jobName is None:
            # Name of executable without path
            jobName = command[:command.find(" ")]
            jobName = jobName[jobName.rfind("/"):]

        print >>f, "#!/bin/bash"
        print >>f, "#   Post this job with `qsub -V $0'"
        print >>f, "#PBS -l nodes=%d:ppn=%d" % (numNodes, numProcsPerNode)
        if wallTime is not None:
            print >>f, "#PBS -l walltime=%d" % wallTime
        print >>f, "#PBS -o %s" % outputDir
        print >>f, "#PBS -N %s" % jobName
        if self.queue is not None:
            print >>f, "#PBS -q %s" % self.queue
        print >>f, "#PBS -j oe"
        print >>f, "#PBS -W umask=02"
        print >>f, "echo \"mpiexec is at: $(which mpiexec)\""
        print >>f, "ulimit -a"
        print >>f, "umask 02"
        print >>f, "echo 'umask: ' $(umask)"
        print >>f, "eups list -s"
        print >>f, "export"
        print >>f, "cd %s" % os.getcwd()
        print >>f, "mpiexec --verbose %s" % command
        f.close()
        os.chmod(script, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
        return script

    def run(command, *args, **kwargs):
        script = self.create(command, *args, **kwargs)
        os.system("qsub -V %s" % script)
        return script


