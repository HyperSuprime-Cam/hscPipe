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
import argparse
from lsst.pex.config import Config
from lsst.pipe.base import ArgumentParser
import hsc.pipe.tasks.overlaps as hscOverlaps

class ParseIntInt(argparse.Action):
    """argparse action callback to parse a '%d,%d' into a tuple"""
    def __call__(self, parser, namespace, values, option_string):
        value = values[-1]
        numbers = value.split(',')
        if len(numbers) != 2:
            parser.error("%s value must be two integers separated by a comma, e.g., '1,2'")
        setattr(namespace, self.dest, tuple(numbers))

def printOverlaps(dataRefList, butler, coadd="deep", tract=None, patch=None, showDataRefs=False, detail=False):
    skyMap = butler.get(namespace.coadd + "Coadd_skyMap")

    tractDictList = [hscOverlaps.getTractPatchList(dataRef, skyMap) for dataRef in dataRefList]
    if showDataRefs:
        print "Tracts for each calexp:"
        for dataRef, tractDict in zip(dataRefList, tractDictList):
            hscOverlaps.printDataRef(dataRef, tractDict, doPatches=detail)
        print

    tractData = hscOverlaps.invert([dataRef.dataId for dataRef in dataRefList], tractDictList)
    if patch is not None:
        if tract is None:
            raise RuntimeError("Tract must be specified if patch is desired.")
        hscOverlaps.printPatch(tract, patch, tractData, doDataId=detail)
    else:
        print "Calexps for each tract/patch:"
        hscOverlaps.printTracts(tractData, tract=tract, doDataId=detail)


if __name__ == "__main__":
    parser = ArgumentParser("hscOverlaps")
    parser.add_id_argument("--id", "raw", help="Data identifier for calexps")
    parser.add_argument("--coadd", type=str, required=True)
    parser.add_argument("--tract", type=int, default=None, help="Tract to print")
    parser.add_argument("--patch", action=ParseIntInt, default=None, help="Patch to print")
    parser.add_argument("--showDataRefs", action="store_true", default=False,
                        help="Show tracts/patches for each dataRef?")
    parser.add_argument("--detail", action="store_true", default=False, help="Print detail?")

    try:
        namespace = parser.parse_args(Config())
    except Exception, e:
        print >> sys.stderr, "Error parsing arguments: %s" % e
        sys.exit(1)

    args = []

    try:
        printOverlaps(namespace.id.refList, namespace.butler, namespace.coadd,
                      namespace.tract, namespace.patch, namespace.showDataRefs, namespace.detail)
    except Exception, e:
        if namespace.doraise:
            raise
        print "Failed: %s" % e
