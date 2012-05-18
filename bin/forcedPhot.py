#!/usr/bin/env python

import argparse, os, sys
import hsc.pipe.tasks.forcedPhot as fp
from hsc.pipe.base import HscArgumentParser

if __name__ == "__main__":
    parser = HscArgumentParser(name="suprimecam")

    try:
        namespace = parser.parse_args(config=fp.HscForcedPhotConfig())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)

    task = fp.HscForcedPhotTask(config=namespace.config)
    butler = namespace.butler
    expIdList = [dataRef.dataId for dataRef in namespace.dataRefList]

    for dataRef in dataRefList:
        if namespace.doRaise:
            task.run(dataRef)
        else:
            try:
                task.run(dataRef)
            except Exception, e:
                task.log.log(task.log.FATAL, "Failed: %s" % e)
