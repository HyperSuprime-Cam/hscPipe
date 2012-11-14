#!/usr/bin/env python

import os
from hsc.pipe.base.pbs import PbsArgumentParser, shCommandFromArgs
from hsc.pipe.tasks.processExposure import ProcessExposureTask


if __name__ == "__main__":
    parent = ProcessExposureTask._makeArgumentParser()
    parser = PbsArgumentParser(description="Reduce frames using PBS.", parent=parent)
    args = parser.parse_args()

    if len(args.parent.dataRefList) == 0:
        print "No frames provided to process"
        exit(1)

    numCcds = sum([len(raft) for raft in args.butler.mapper.camera])
    command = "python %s/bin/processExposure.py %s" % (os.environ['HSCPIPE_DIR'],
                                                       shCommandFromArgs(args.leftover))
    args.pbs.run(command, repeats=len(args.frame), threads=numCcds)
