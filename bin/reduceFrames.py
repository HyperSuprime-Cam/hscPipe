#!/usr/bin/env python

import os
import lsst.afw.cameraGeom as cameraGeom
from hsc.pipe.base.pbs import PbsArgumentParser, shCommandFromArgs
from hsc.pipe.tasks.processExposure import ProcessExposureTask


if __name__ == "__main__":
    processParser = ProcessExposureTask._makeArgumentParser(add_help=False)
    pbsParser = PbsArgumentParser(description="Reduce frames using PBS.", parent=processParser)
    args = pbsParser.parse_args(config=ProcessExposureTask.ConfigClass())

    numExps = len(args.parent.dataRefList)
    if numExps == 0:
        print "No frames provided to process"
        exit(1)

    numCcds = sum([sum([1 for ccd in cameraGeom.cast_Raft(raft)])
                   for raft in args.parent.butler.mapper.camera])
    command = "python %s/bin/processExposure.py %s" % (os.environ['HSCPIPE_DIR'],
                                                       shCommandFromArgs(args.leftover))
    args.pbs.run(command, repeats=numExps, threads=numCcds)
