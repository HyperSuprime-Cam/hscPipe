#!/usr/bin/env python

import argparse, os, sys
import hsc.pipe.tasks.forcedPhot as fp
from lsst.pipe.base import ArgumentParser

class OutputAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if namespace.rerun:
            raise argparse.ArgumentTypeError("Please specify --output or --rerun, but not both")

        namespace.outPath = values

class RerunAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        """We can't just parse the arguments and reset namespace.outPath as the Mapper's been
        instantiated before we get a chance"""

        if namespace.outPath:
            raise argparse.ArgumentTypeError("Please specify --output or --rerun, but not both")

        envar = "SUPRIME_DATA_DIR"
        if os.environ.has_key(envar):
            namespace.rerun = values
            namespace.outPath = os.path.join(os.environ[envar], "SUPA", "rerun", namespace.rerun)
            if not os.path.exists(namespace.outPath):
                os.makedirs(namespace.outPath) # should be in butler
        else:
            raise argparse.ArgumentTypeError("You must define $%s to use --rerun XXX" % envar)

if __name__ == "__main__":
    try:
        parser = ArgumentParser(name="suprimecam", conflict_handler='resolve') # new style
    except TypeError:
        parser = ArgumentParser(conflict_handler='resolve') # old style

    parser.add_argument('--output', type=str, dest="outPath", default=None, help="output root directory",
                        action=OutputAction)
    parser.add_argument('--rerun', type=str, default=None, help='Desired rerun (overrides --output)',
                        action=RerunAction)

    parser.add_argument('--stack', type=int, optional=False, help='Stack identifier')
    parser.add_argument('--filter', type=str, optional=False, help='Filter name')

    try:
        namespace = parser.parse_args(config=fp.ForcedPhotConfig())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)
            
    task = fp.ForcedPhotTask(config=namespace.config)
    stackId = {'stack': namespace.stack, 'filter': namespace.filter}
    expIdList = [dataRef.dataId for dataRef in namespace.dataRefList]

    if namespace.doRaise:
        task.run(butler, stackId, expIdList)
    else:
        try:
            task.run(butler, stackId, expIdList)
        except Exception, e:
            task.log.log(task.log.FATAL, "Failed on dataId=%s: %s" % (sensorRef.dataId, e))
