#!/usr/bin/env python

import os
import hsc.pipe.base.pbs as hscPbs


if __name__ == "__main__":
    parser = hscPbs.PbsArgumentParser(usage="""
    Reduce frames using PBS.

    Command line:
        reduceFrames [OPTIONS] FRAME_ID [FRAME_ID ...]

    where FRAME_ID are integers identifying frames.
    """)
    parser.add_argument("-i", "--instrument", dest="instrument", help="Instrument name", default="hsc")
    parser.add_argument("frame", nargs='*', help="Frame numbers to reduce")
    pbs, args = parser.parse_args()

    if len(args.frame) == 0:
        print "No frames provided to process"
        exit(1)

    command = "python %s/bin/processExposure.py %s %s" % (os.environ['HSCPIPE_DIR'], args.instrument,
                                                          " ".join(args.frame))
    pbs.run(command, repeats=len(args.frame))
