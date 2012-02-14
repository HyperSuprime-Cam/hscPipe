#!/usr/bin/env python

import os
import os.path
import stat
import sys
import tempfile
import argparse


class PbsArgumentParser(argparse.ArgumentParser):
    """An argument parser to get relevant parameters for PBS."""
    def __init__(self, *args, **kwargs):
        super(PbsArgumentParser, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--queue", dest="queue", help="PBS queue name")
        self.add_argument("-r", "--rerun", dest="rerun", help="Rerun name")
        self.add_argument("-j", "--job", dest="job", help="Job name")
        self.add_argument("-n", "--nodes", dest="nodes", help="Number of nodes")
        self.add_argument("-p", "--procs", dest="procs", help="Number of processors per node")
        self.add_argument("-t", "--time", dest="time", help="Expected execution time per processor (sec)")
        self.add_argument("-o", "--output", dest="output", help="Output directory")
        self.add_argument("-N", "--dry-run", dest="dryrun", default=False, action="store_true",
                          help="Dry run?")

    def getPbs(self):
        return Pbs(outputDir=self.output, numNodes=self.nodes, numProcsPerNode=self.procs,
                   queue=self.queue, jobName=self.job, wallTime=self.time, dryrun=self.dryrun)


class Pbs(object):
    def __init__(self, outputDir=None, numNodes=1, numProcsPerNode=1, queue=None, jobName=None, time=None,
                 dryrun=False):
        self.outputDir = outputDir
        self.numNodes = numNodes
        self.numProcsPerNode = numProcsPerNode
        self.queue = queue
        self.jobName = jobName
        self.time = wallTime

    def create(command, repeats=1, time=None, numNodes=None, numProcsPerNode=None, jobName=None):
        if time is None:
            time = self.time
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
        if time is not None:
            wallTime = repeats * time * numNodes * numProcsPerNode
            print >>f, "#PBS -l walltime=%d" % wallTime
        if self.outputDir is not None:
            print >>f, "#PBS -o %s" % self.outputDir
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
        command = "qsub -V %s" % script
        if self.dryrun:
            print "Would run: %s" % command
        else:
            os.system(command)
        return script


