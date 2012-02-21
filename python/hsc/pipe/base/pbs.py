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
        self.add_argument("-j", "--job", dest="job", help="Job name")
        self.add_argument("-n", "--nodes", dest="nodes", type=int, help="Number of nodes", required=True)
        self.add_argument("-p", "--procs", dest="procs", type=int, help="Number of processors per node",
                          required=True)
        self.add_argument("-t", "--time", dest="time", type=float,
                          help="Expected execution time per processor (sec)")
        self.add_argument("-o", "--output", dest="output", help="Output directory")
        self.add_argument("-N", "--dry-run", dest="dryrun", default=False, action="store_true",
                          help="Dry run?")

    def parse_args(self, *args, **kwargs):
        args = super(PbsArgumentParser, self).parse_args(*args, **kwargs)
        pbs = Pbs(outputDir=args.output, numNodes=args.nodes, numProcsPerNode=args.procs,
                  queue=args.queue, jobName=args.job, time=args.time, dryrun=args.dryrun)
        return pbs, args

class Pbs(object):
    def __init__(self, outputDir=None, numNodes=1, numProcsPerNode=1, queue=None, jobName=None, time=None,
                 dryrun=False):
        self.outputDir = outputDir
        self.numNodes = numNodes
        self.numProcsPerNode = numProcsPerNode
        self.queue = queue
        self.jobName = jobName
        self.time = time
        self.dryrun = dryrun

    def create(self, command, repeats=1, time=None, numNodes=None, numProcsPerNode=None, jobName=None,
               threads=None):
        if time is None:
            time = self.time
        if numNodes is None:
            numNodes = self.numNodes
        if numProcsPerNode is None:
            numProcsPerNode = self.numProcsPerNode
        if jobName is None:
            jobName = self.jobName
        if threads is None:
            threads = numNodes * numProcsPerNode
        threads = min(threads, numNodes * numProcsPerNode)

        fd, script = tempfile.mkstemp()
        f = os.fdopen(fd, "w")

        if numNodes is None or numProcsPerNode is None:
            raise RuntimeError("numNodes (%s) or numProcsPerNode (%s) is not specified" %
                               (numNodes, numProcsPerNode))

        assert numNodes is not None and numProcsPerNode is not None
        if jobName is None:
            # Name of executable without path
            jobName = command[:command.find(" ")]
            jobName = jobName[jobName.rfind("/"):]

        print >>f, "#!/bin/bash"
        print >>f, "#   Post this job with `qsub -V $0'"
        print >>f, "#PBS -l nodes=%d:ppn=%d" % (numNodes, numProcsPerNode)
        if time is not None:
            wallTime = repeats * time / threads
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

    def run(self, command, *args, **kwargs):
        script = self.create(command, *args, **kwargs)
        command = "qsub -V %s" % script
        if self.dryrun:
            print "Would run: %s" % command
        else:
            os.system(command)
        return script


