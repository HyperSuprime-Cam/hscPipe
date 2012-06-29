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
import hsc.pipe.tasks.monkeypatch
from lsst.pipe.base import ArgumentParser
from hsc.meas.mosaic.task import HscCoaddTask as TaskClass

if __name__ == "__main__":
    coaddName = "deep"
    parser = ArgumentParser("hscCoadd", datasetType="calexp")
    parser.add_argument("--filter", type=str, required=True)
    parser.add_argument("--tract", type=int, required=True)
    parser.add_argument("--patch", type=str, required=True)

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)

    task = TaskClass(config=namespace.config)
    coaddId = {'filter': namespace.filter,
               'tract': namespace.tract,
               'patch': namespace.patch,
               }
    coaddRefList = list(namespace.butler.subset(datasetType=namespace.config.coaddName + "Coadd_skyMap",
                        level=None, dataId=coaddId))
    if len(coaddRefList) != 1:
        raise RuntimeError("Non-specific coadd reference list: %s" % coaddRefList)

    
    if namespace.doRaise:
        task.run(coaddRefList[0], namespace.dataRefList)
    else:
        try:
            task.run(coaddRefList[0], namespace.dataRefList)
        except Exception, e:
            task.log.log(task.log.FATAL, "Failed: %s" % e)
