#!/usr/bin/env python

import os
import hsc.pipe.base.pbs as hscPbs

if __name__ == "__main__":
    parser = hscPbs.PbsArgumentParser(usage="""
    Stack frames using PBS.

    Command line:
        stackFrames [OPTIONS] FIELD_NAME FILTER_NAME
    """)
    parser.add_argument("-i", "--instrument", dest="instrument", help="Instrument name", default="hsc")
    parser.add_argument("-r", "--rerun", dest="rerun", help="Rerun name", default=os.environ["LOGNAME"])
    parser.add_argument("--wcs", dest="wcs", help="WCS for stack")
    parser.add_argument("-m", "--doMatchPsf", default=False, action='store_true',
                        help="match PSFs before stacking?")
    parser.add_argument("field", help="Field name in butler for stack")
    parser.add_argument("filter", help="Filter name for stack")
    pbs, args = parser.parse_args()
    assert args.field and args.filter, "Field and filter not provided"

    command = "python %s/bin/stackExposures.py" % os.environ["HSCPIPE_DIR"]
    command += " --instrument=%s --rerun=%s" % (args.instrument, args.rerun)
    if args.output == None:
        command += " --program=%s --filter=%s --workDirRoot=%s" % (args.field, args.filter, ".")
    else:
        command += " --program=%s --filter=%s --workDirRoot=%s" % (args.field, args.filter, args.output)
    if args.doMatchPsf:
        command += " --doMatchPsf"
    if args.wcs is not None:
        command += " --destWcs=%s" % args.wcs
    command += " --pid=%d" % os.getpid()

    pbs.run(command)
