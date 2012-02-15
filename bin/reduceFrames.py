#!/usr/bin/env python

import hsc.pipe.base.pbs as hscPbs


if __name__ == "__main__":
    parser = hscPbs.PbsArgumentParser(usage="""
    Reduce frames using PBS.

    Command line:
        reduceFrames [OPTIONS] FRAME_ID [FRAME_ID ...]

    where FRAME_ID are integers identifying frames.
    """)
    parser.add_argument("-i", "--instrument", dest="instrument", help="Instrument name", default="hsc")
    parser.add_argument("frame", nargs='*', type=int, help="Frame numbers to reduce")
    pbs, args = parser.parse_args()

    command = "python %s/bin/processExposure.py %s %s" % (os.environ['HSCPIPE_DIR'], args.instrument,
                                                          " ".join(args.frame))
    pbs.run(command, repeats=len(args.frame))
