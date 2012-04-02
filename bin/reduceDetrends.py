#!/usr/bin/env python

import os
import hsc.pipe.base.pbs as hscPbs
import hsc.pipe.base.camera as hscCamera


if __name__ == "__main__":
    parser = hscPbs.PbsArgumentParser(usage="""
    Reduce detrend frames using PBS.

    Command line:
        reduceDetrends [OPTIONS] FRAME_ID [FRAME_ID ...]

    where FRAME_ID are integers identifying frames.
    """)
    parser.add_argument("-i", "--instrument", dest="instrument", help="Instrument name", default="hsc")
    parser.add_argument("-r", "--rerun", dest="rerun", help="Rerun name", default=os.getlogin())
    parser.add_argument("-d", "--detrend", type=str, choices=["bias", "dark", "flat", "fringe", "mask"],
                        help="Detrend type")
    parser.add_argument("-o", "--out", type=str, help="Pattern for output name, with %%d for ccd number")
    parser.add_argument("frame", nargs='*', help="Frame numbers to reduce")
    pbs, args = parser.parse_args()

    if len(args.frame) == 0:
        print "No frames provided to process"
        exit(1)

    command = "python %s/bin/processDetrends.py --instrument %s --rerun %s --detrend %s --out %s" % \
        (os.environ['HSCPIPE_DIR'], args.instrument, args.rerun, args.detrend, args.out " ".join(args.frame))
    pbs.run(command, repeats=len(args.frame)*hscCamera.getNumCcds(args.instrument))
