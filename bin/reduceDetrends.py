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
    parser.add_argument("-r", "--rerun", dest="rerun", help="Rerun name", default=os.environ["LOGNAME"])
    parser.add_argument("-d", "--detrend", type=str, choices=["bias", "dark", "flat", "fringe", "mask"],
                        help="Detrend type", required=True)
    parser.add_argument("-O", "--out", type=str, help="Pattern for output name, with %%d for ccd number")
    parser.add_argument("frame", nargs='*', help="Frame numbers to reduce")
    pbs, args = parser.parse_args()

    if len(args.frame) == 0:
        print "No frames provided to process"
        exit(1)

    command = "python %s/bin/processDetrends.py" % os.environ['HSCPIPE_DIR']
    command += " --instrument %s" % args.instrument
    command += " --rerun %s" % args.rerun
    command += " --detrend %s" % args.detrend
    if args.out is not None:
        command += " --out %s" % args.out
    command += " " + " ".join(args.frame)
    pbs.run(command, repeats=len(args.frame)*hscCamera.getNumCcds(args.instrument))
