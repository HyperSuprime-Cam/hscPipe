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

import os
import numpy

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import hsc.pipe.base as hscPipeBase
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
from hsc.pipe.tasks.processStack import StackArgumentParser

class MergeStackSourcesConfig(pexConfig.Config):
    patch = pexConfig.Field(dtype=int, default=999999, optional=False, doc="patch dataId key for output")

class MergeStackSourcesTask(pipeBase.CmdLineTask):
    """Process a Stack image
    
    """
    ConfigClass = MergeStackSourcesConfig
    _DefaultName = "mergeStackSources"

    # Defining these to return None turns off warnings about persisting configs and metadata
    def _getConfigName(self): return None
    def _getMetadataName(self): return None

    def __init__(self, **kwargs):
        pipeBase.CmdLineTask.__init__(self, **kwargs)
        self.outputCat = None
        self.outputWcs = None

    def loadWcs(self, dataRef, dataId=None):
        butler = dataRef.butlerSubset.butler
        if dataId is None:
            dataId = dataRef.dataId
        try:
            exp = butler.get("stack_calexp", dataId, immediate=True)
        except:
            exp = butler.get("stack", dataId, immediate=True)
        return exp.getWcs()

    @pipeBase.timeMethod
    def run(self, dataRef):
        if self.outputWcs is None:
            self.log.info("Loading output Wcs")
            self.outputId = dataRef.dataId.copy()
            self.outputId['patch'] = self.config.patch
            self.outputWcs = self.loadWcs(dataRef, dataId=self.outputId)
        self.log.info("Loading sources for %s" % dataRef.dataId)
        inputCat = dataRef.get("stack_src", immediate=True)
        inputWcs = self.loadWcs(dataRef)
        self.log.info("Shifting sources for %s" % dataRef.dataId)
        offsetD = self.outputWcs.getPixelOrigin() - inputWcs.getPixelOrigin()
        offsetI = afwGeom.Extent2I(int(offsetD.getX()), int(offsetD.getY()))
        for record in inputCat:
            footprint = record.getFootprint()
            footprint.shift(offsetI)
            record.setFootprint(footprint)
            # Should set other centroids too, but Jim is in a hurry and he knows we don't
            # care about them for forced photometry
            record.set(inputCat.table.getCentroidKey(), record.getCentroid() + offsetD)
        if self.outputCat is None:
            self.outputCat = afwTable.SourceCatalog(inputCat.table.clone())
        self.log.info("Transferring sources for %s" % dataRef.dataId)
        self.outputCat.extend(inputCat, deep=True)

    @pipeBase.timeMethod
    def finish(self, butler):
        butler.put(self.outputCat, "stack_src", self.outputId)

    @classmethod
    def _makeArgumentParser(cls):
        return StackArgumentParser(name=cls._DefaultName, datasetType="stack", defaultToUserName=False)
