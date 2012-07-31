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
import sys
from hsc.pipe.base import HscArgumentParser
from hsc.meas.mosaic.task import MosaicTask as TaskClass

if __name__ == "__main__":
    coaddName = "deep"
    parser = HscArgumentParser(datasetType="calexp")
    parser.add_argument("--filter", type=str)
    parser.add_argument("--tract", type=int)

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)

    task = TaskClass(config=namespace.config)
    butler = namespace.butler
    tractId = {}
    if namespace.tract is not None:
        tractId['tract'] = namespace.tract
    if namespace.doRaise:
        task.run(butler, namespace.dataRefList, tractId, coaddName)
    else:
        try:
            task.run(butler, namespace.dataRefList, tractId, coaddName)
        except Exception, e:
            task.log.log(task.log.FATAL, "Failed on dataId=%s: %s" % (dataId, e))