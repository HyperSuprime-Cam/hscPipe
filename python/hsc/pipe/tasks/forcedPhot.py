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

from lsst.pipe.tasks.forcedPhot import ForcedPhotTask, ForcedPhotConfig, ReferencesTask, ReferencesConfig
from lsst.pex.config import Field, ConfigurableField

class SubaruReferencesConfig(ReferencesConfig):
    filter = Field(dtype=str, doc="Filter name of stack for references", optional=False)

class SubaruReferencesTask(ReferencesTask):
    ConfigClass = SubaruReferencesConfig

    def getReferenceSources(self, dataRef, exposure):
        """Get Wcs for reference sources"""
        butler = dataRef.butlerSubset.butler
        dataId = dataRef.dataId
        filterName = self.config.filter if self.config.filter is not None else dataId['filter']
        stackId = {'stack': dataId['pointing'],
                   'patch': 999999,
                   'filter': filterName,
                   }
        return butler.get("stack_src", stackId, immediate=True)

    def getReferenceWcs(self, dataRef, exposure):
        """Get reference sources on (or close to) exposure"""
        butler = dataRef.butlerSubset.butler
        dataId = dataRef.dataId
        filterName = self.config.filter if self.config.filter is not None else dataId['filter']
        stackId = {'stack': dataId['pointing'],
                   'patch': 999999,
                   'filter': filterName,
                   }
        calexp = butler.get("stack_calexp", stackId, immediate=True)
        return calexp.getWcs()

class SubaruForcedPhotTask(ForcedPhotTask):

    def readInputs(self, dataRef, *args, **kwargs):
        struct = super(SubaruForcedPhotTask, self).readInputs(dataRef, *args, **kwargs)
        try:
            wcsExp = dataRef.get("wcs")
            struct.exposure.setWcs(wcsExp.getWcs())
        except Exception, e:
            self.log.log(self.log.WARN, "No WCS tweak available (%s); not updating WCS." % e)

        return struct
