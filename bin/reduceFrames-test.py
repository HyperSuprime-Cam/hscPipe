#!/usr/bin/env python

import os
import hsc.pipe.base.parallel as hscParallel


if __name__ == "__main__":
    parser = hscParallel.BatchArgumentParser(usage="""
    Reduce frames using batch system.

    Command line:
        reduceFrames [OPTIONS] FRAME_ID [FRAME_ID ...]

    where FRAME_ID are integers identifying frames.
    """)
    parser.add_argument("-i", "--instrument", dest="instrument", help="Instrument name", default="hsc")
    parser.add_argument("-r", "--rerun", dest="rerun", help="Rerun name", default=os.environ["LOGNAME"])
    parser.add_argument("-I", "--input-data-dir", dest="root", help="directory of input data (e.g., /data/Subaru/HSC)")
    parser.add_argument("-O", "--output-data-dir", dest="outputRoot", help="directory to output data (e.g., /data/rerun/XXX)")
    parser.add_argument("-C", "--calib-dir", dest="calibRoot", help="directory of calibration data (e.g., /data/Subaru/HSC/CALIB)")
    parser.add_argument("-c", "--configfile", dest="configFile", help="config override file")
    parser.add_argument("frame", nargs='*', help="Frame numbers to reduce")
    batch, args = parser.parse_args()

    if len(args.frame) == 0:
        print "No frames provided to process"
        exit(1)

    command = "python %s/bin/processExposure-test.py %s %s %s %s %s %s %s" % (os.environ['HSCPIPE_DIR'], args.instrument,
                                                                   args.rerun,
                                                                   args.root, args.outputRoot,args.calibRoot, 
                                                                      args.configFile,
                                                                   " ".join(args.frame))
    batch.run(command, repeats=len(args.frame), threads=10 if args.instrument == "suprimecam" else None)
