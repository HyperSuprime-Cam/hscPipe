#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import argparse, os, sys
from lsst.pipe.base import ArgumentParser
from hsc.pipe.tasks.processCcd import SuprimeCamProcessCcdTask as TaskClass

# FH added for QA output
# import hsc.onsite.qa.fitsthumb as QaFitsthumb
# import hsc.onsite.qa.measSeeingQa as QaSeeing


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

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)
            
    task = TaskClass(config=namespace.config)
    for sensorRef in namespace.dataRefList:
        if namespace.doRaise:
            task.run(sensorRef)
        else:
            try:
                task.run(sensorRef)
            except Exception, e:
                task.log.log(task.log.FATAL, "Failed on dataId=%s: %s" % (sensorRef.dataId, e))
