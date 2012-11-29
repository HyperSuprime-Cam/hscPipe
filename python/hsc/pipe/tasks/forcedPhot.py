#!/usr/bin/env python

import lsst.daf.base as dafBase
import lsst.daf.persistence as dafPersist
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.detection as afwDet
import lsst.ip.isr as ipIsr
import lsst.pipe.tasks.calibrate as ptCal
import lsst.meas.algorithms as measAlg
import lsst.meas.mosaic as measMosaic

from lsst.pipe.tasks.forcedPhot import ForcedPhotTask, ForcedPhotConfig, ReferencesTask, ReferencesConfig
from lsst.pex.config import Field, ConfigurableField

class SubaruReferencesConfig(ReferencesConfig):
    filter = Field(dtype=str, doc="Filter name of stack for references", optional=False)

class SubaruReferencesTask(ReferencesTask):
    ConfigClass = SubaruReferencesConfig

    def getReferenceSources(self, dataRef, exposure):
        """Get reference sources"""
        if not hasattr(self, "_referenceSources"):
            self._referenceSources = {}
        butler = dataRef.butlerSubset.butler
        dataId = dataRef.dataId
        filterName = self.config.filter if self.config.filter is not None else dataId['filter']
        key = (dataId["pointing"], filterName)
        if key not in self._referenceSources:
            stackId = {'stack': dataId['pointing'],
                       'patch': 999999,
                       'filter': filterName,
                       }
            self._referenceSources[key] = butler.get("stack_src", stackId, immediate=True)
        return self._referenceSources[key]

    def getReferenceWcs(self, dataRef, exposure):
        """Get Wcs for reference sources"""
        if not hasattr(self, "_referenceWcs"):
            self._referenceWcs = {}
        butler = dataRef.butlerSubset.butler
        dataId = dataRef.dataId
        filterName = self.config.filter if self.config.filter is not None else dataId['filter']
        key = (dataId["pointing"], filterName)
        if key not in self._referenceWcs:
            stackId = {'stack': dataId['pointing'],
                       'patch': 999999,
                       'filter': filterName,
                       }
            calexp = butler.get("stack_calexp", stackId, immediate=True)
            self._referenceWcs[key] = calexp.getWcs()
        return self._referenceWcs[key]

class SubaruForcedPhotConfig(ForcedPhotTask.ConfigClass):
    useMosaicWcs = Field(dtype=bool, optional=False, default=True,
                         doc="Whether to use the Wcs saved by meas_mosaic.")
    useMosaicFluxFit = Field(dtype=bool, optional=False, default=True,
                             doc="Whether to apply the FluxFitParams saved by meas_mosaic.")

class SubaruForcedPhotTask(ForcedPhotTask):

    ConfigClass = SubaruForcedPhotConfig

    def readInputs(self, dataRef, *args, **kwargs):
        struct = super(SubaruForcedPhotTask, self).readInputs(dataRef, *args, **kwargs)
        if self.config.useMosaicWcs:
            wcsExp = dataRef.get("wcs", immediate=True)
            struct.exposure.setWcs(wcsExp.getWcs())
        if self.config.useMosaicFluxFit:
            ffp = measMosaic.FluxFitParams(dataRef.get("fcr_md", immediate=True))
            mi = struct.exposure.getMaskedImage()
            fcor = measMosaic.getFCorImg(ffp, mi.getWidth(), mi.getHeight())
            mi *= fcor
        return struct
