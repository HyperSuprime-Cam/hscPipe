#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.calibrate as ptCalibrate
import lsst.meas.algorithms as measAlg
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.photocal as photocal
from lsst.pipe.tasks.repair import RepairTask
from lsst.pipe.tasks.measurePsf import MeasurePsfTask
import hsc.pipe.tasks.astrometry as hscAstrom


class SubaruCalibrateConfig(ptCalibrate.CalibrateConfig):
    astrometry = pexConfig.ConfigField(dtype = hscAstrom.SubaruAstrometryConfig, doc = "HSC calibration")

class SubaruCalibrateTask(ptCalibrate.CalibrateTask):
    ConfigClass = SubaruCalibrateConfig

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        self.makeSubtask("repair", RepairTask)
        self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        self.makeSubtask("initialMeasurement", measAlg.SourceMeasurementTask,
                         schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("measurePsf", MeasurePsfTask, schema=self.schema)
        self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                         schema=self.schema, algMetadata=self.algMetadata)
        self.makeSubtask("astrometry", hscAstrom.SubaruAstrometryTask, schema=self.schema)
        self.makeSubtask("photocal", photocal.PhotoCalTask, schema=self.schema)

    def run(self, exposure, *args, **kwargs):
        results = super(SubaruCalibrateTask, self).run(exposure, *args, **kwargs)

        photocal = results.photocal
        magZero = photocal.zp - 2.5 * math.log10(exposure.getCalib().getExptime()) # convert to (mag/sec/adu)
        self.metadata.set('MAGZERO', magZero)
        self.metadata.set('MAGZERO_RMS', photocal.sigma)
        self.metadata.set('MAGZERO_NOBJ', photocal.ngood)
        self.metadata.set('COLORTERM1', 0.0)
        self.metadata.set('COLORTERM2', 0.0)
        self.metadata.set('COLORTERM3', 0.0)

        return results
