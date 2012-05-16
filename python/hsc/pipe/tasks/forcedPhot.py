#!/usr/bin/env python

import lsst.daf.persistence as dafPersist
import lsst.afw.table as afwTable
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.pipe.tasks.calibrate as ptCal
import lsst.meas.algorithms as measAlg

from lsst.pipe.tasks.forcedPhot import ForcedPhotTask, ReferencesTask, ReferencesConfig
from lsst.pex.config import ConfigurableField

class HscReferencesConfig(ReferencesConfig):
    calibrate = ConfigurableField(target=ptCal.CalibrateConfig, doc="Configuration for calibration of stack")
    detection = ConfigurableField(target=measAlg.SourceDetectionConfig,
                                  doc="Configuration for detection on stack")

class HscReferencesTask(ReferencesTask):
    def __init__(self, *args, **kwargs):
        super(ReferencesTask, self).__init__(*args, **kwargs)
        self.makeSubtask("calibrate")
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("detection", schema=self.schema)

    def getReferences(self, dataRef, exposure):
        """Get reference sources on (or close to) exposure"""
        butler = dataRef.butler
        dataId = dataRef.dataId
        stackId = {'stack': dataId['pointing'],
                   'patch': 999999,
                   'filter': dataId['filter'],
                   }
        try:
            sources = foilReadProxy(butler.get("stacksources", stackId))
        except Exception, e:
            self.log.log(self.log.INFO, "Stack products not available (%s); attempting to create" % e)
            sources = measureStack(butler, stackId)
        
        return measureStack(butler, stackId)

    def measureStack(self, butler, dataId):
        exposure = butler.get("stack", dataId)
        self.interpolateNans(exposure, afwDet.createPsf("DoubleGaussian", 15, 15, 1.0))       
        calib = self.calibrate.run(exposure)
        exposure = calib.exposure

        det = self.detect(exposure)
        sources = det.sources
        self.measurement.run(exposure, sources, apCorr=calib.apCorr)
        butler.put(sources, "stacksources", dataId)
        butler.put(exposure.getPsf(), "stackpsf", dataId)

        return sources


    def interpolateNans(self, exposure, psf):
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN")
        nanMasker = ipIsr.UnmaskedNanCounterF()
        nanMasker.apply(exposure.getMaskedImage())
        nans = ipIsr.Isr().getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN')
        self.log.log(self.log.INFO, "Interpolating over %d NANs" % len(nans))
        measAlg.interpolateOverDefects(exposure.getMaskedImage(), psf, nans, 0.0)
        
    def detect(self, exposure):
        table = afwTable.SourceTable.make(self.schema)
        table.setMetadata(self.algMetadata)
        return self.detection.makeSourceCatalog(table, exposure)


class HscForcedPhotTask(ForcedPhotTask):
    def readInputs(self, dataRef, *args, **kwargs):
        struct = super(HscForcedPhotTask, self).readInputs(dataRef, *args, **kwargs)
        struct.exposure.setWcs(dataRef.get("wcs"))
        return struct


def foilReadProxy(obj):
    if isinstance(obj, dafPersist.ReadProxy):
        obj = obj.__subject__
    return obj
