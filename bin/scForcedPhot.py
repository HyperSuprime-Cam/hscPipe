#!/usr/bin/env python

import argparse, os, sys
import hsc.pipe.tasks.forcedPhot as fp
from hsc.pipe.base import HscArgumentParser

if __name__ == "__main__":
    parser = HscArgumentParser()

    try:
        namespace = parser.parse_args(config=fp.HscForcedPhotConfig())
    except Exception, e:
        if "--doraise" in sys.argv:
            raise
        print >> sys.stderr, e
        sys.exit(1)

    task = fp.HscForcedPhotTask(config=namespace.config)

    for dataRef in namespace.dataRefList:
        if namespace.doRaise:
            task.run(dataRef)
        else:
            try:
                task.run(dataRef)
            except Exception, e:
                task.log.log(task.log.FATAL, "Failed: %s" % e)
