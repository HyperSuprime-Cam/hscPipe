#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.calibrate as ptCalibrate
import lsst.meas.algorithms as measAlg
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
###import lsst.meas.photocal as photocal
from lsst.pipe.tasks.repair import RepairTask
from lsst.pipe.tasks.measurePsf import MeasurePsfTask
import hsc.pipe.tasks.astrometry as hscAstrom

import hsc.pipe.tasks.photocalQa as photocalQa
import math

class HscCalibrateConfig(ptCalibrate.CalibrateConfig):
    astrometry = pexConfig.ConfigField(dtype = hscAstrom.HscAstrometryConfig, doc = "HSC calibration")

class HscCalibrateTask(ptCalibrate.CalibrateTask):
    ConfigClass = HscCalibrateConfig

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
        self.makeSubtask("astrometry", hscAstrom.HscAstrometryTask, schema=self.schema)
##== FH modified for QA output
##        self.makeSubtask("photocal", photocal.PhotoCalTask, schema=self.schema)
        self.makeSubtask("photocal", photocalQa.PhotoCalTask, schema=self.schema)        

##== FH modified for QA output
    @pipeBase.timeMethod
    def run(self, exposure, defects=None):
        """Calibrate an exposure: measure PSF, subtract background, measure astrometry and photometry

        @param[in,out]  exposure   Exposure to calibrate; measured PSF will be installed there as well
        @param[in]      defects    List of defects on exposure
        @return a pipeBase.Struct with fields:
        - psf: Point spread function
        - apCorr: Aperture correction
        - sources: Sources used in calibration
        - matches: Astrometric matches
        - matchMeta: Metadata for astrometric matches
        """
        assert exposure is not None, "No exposure provided"

        self.installInitialPsf(exposure)

        keepCRs = True                  # At least until we know the PSF
        self.repair.run(exposure, defects=defects, keepCRs=keepCRs)
        self.display('repair', exposure=exposure)

        if self.config.doBackground:
            with self.timer("background"):
                bg, exposure = measAlg.estimateBackground(exposure, self.config.background, subtract=True)
                del bg

            self.display('background', exposure=exposure)
        
        table = afwTable.SourceTable.make(self.schema) # TODO: custom IdFactory for globally unique IDs
        table.setMetadata(self.algMetadata)
        detRet = self.detection.makeSourceCatalog(table, exposure)
        sources = detRet.sources

        if self.config.doPsf:
            self.initialMeasurement.measure(exposure, sources)
            psfRet = self.measurePsf.run(exposure, sources)
            cellSet = psfRet.cellSet
            psf = psfRet.psf
        else:
            psf, cellSet = None, None

        # Wash, rinse, repeat with proper PSF

        if self.config.doPsf:
            self.repair.run(exposure, defects=defects, keepCRs=None)
            self.display('repair', exposure=exposure)

        if self.config.doBackground:   # is repeating this necessary?  (does background depend on PSF model?)
            with self.timer("background"):
                # Subtract background
                background, exposure = measAlg.estimateBackground(
                    exposure, self.config.background, subtract=True)
                self.log.log(self.log.INFO, "Fit and subtracted background")

            self.display('background', exposure=exposure)

        if self.config.doComputeApCorr or self.config.doAstrometry or self.config.doPhotoCal:
            self.measurement.measure(exposure, sources)   # don't use run, because we don't have apCorr yet

        if self.config.doComputeApCorr:
            assert(self.config.doPsf)
            apCorr = self.computeApCorr(exposure, cellSet)
        else:
            apCorr = None

        if self.measurement.config.doApplyApCorr:
            assert(apCorr is not None)
            self.measurement.applyApCorr(sources, apCorr)

        if self.config.doAstrometry:
            astromRet = self.astrometry.run(exposure, sources)
            matches = astromRet.matches
            matchMeta = astromRet.matchMeta
        else:
            matches, matchMeta = None, None

        if self.config.doPhotoCal:
            assert(matches is not None)
            photocalRet = self.photocal.run(matches)
            zp = photocalRet.photocal
            self.log.log(self.log.INFO, "Photometric zero-point: %f" % zp.getMag(1.0))
            exposure.getCalib().setFluxMag0(zp.getFlux(0))

            ##== FH: added for QA output
            magZero = photocalRet.zp
            magSigma = photocalRet.sigma
            nRef = photocalRet.ngood
            calib = exposure.getCalib()
            calib.setFluxMag0(zp.getFlux(0))
            exptime = calib.getExptime()
            self.log.log(self.log.INFO, "QA photocal: exptime: %f (sec)" % exptime)
            self.log.log(self.log.INFO, "QA photocal: magzero: %f" % magZero)
            self.log.log(self.log.INFO, "QA photocal: magzeroRms: %f" % magSigma)
            self.log.log(self.log.INFO, "QA photocal: magzeroNref: %d" % nRef)                        
            
            magZero = magZero - 2.5 * math.log10(exptime) # convert to (mag/sec/adu)
            metadata = exposure.getMetadata()
            metadata.set('MAGZERO', magZero)
            metadata.set('MAGZERO_RMS', magSigma)
            metadata.set('MAGZERO_NOBJ', nRef)
            metadata.set('COLORTERM1', 0.0)
            metadata.set('COLORTERM2', 0.0)
            metadata.set('COLORTERM3', 0.0)

        self.display('calibrate', exposure=exposure, sources=sources, matches=matches)

        return pipeBase.Struct(
            exposure = exposure,
            psf = psf,
            apCorr = apCorr,
            sources = sources,
            matches = matches,
            matchMeta = matchMeta,
        )

