#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDet
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
from lsst.ip.isr import IsrTask
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.processCcd as ptProcessCcd
import hsc.pipe.tasks.astrometry as hscAstrom
import hsc.pipe.tasks.suprimecam as hscSuprimeCam
import hsc.pipe.tasks.qaHscSuprimeCamTask as qaHscSuprimeCamIsrTask
import hsc.pipe.tasks.calibrate as hscCalibrate
import hsc.pipe.tasks.hscDc2 as hscDc2

##== FH added for QA output
import hsc.onsite.qa.measSeeingQa as QaSeeing

#== FH changed for QA output
class QaFlatness(pexConfig.Config):
    meshX = pexConfig.Field(
        dtype = int,
        doc = 'Mesh size in X (pix) to calculate count statistics',
        default = 256,
        )
    meshY = pexConfig.Field(
        dtype = int,
        doc = 'Mesh size in Y (pix) to calculate count statistics',
        default = 256,
        )
    doClip = pexConfig.Field(
        dtype = bool,
        doc = 'Do we clip outliers in calculate count statistics?',
        default = True,
        )
    clipSigma = pexConfig.Field(
        dtype = float,
        doc = 'How many sigma is used to clip outliers in calculate count statistics?',
        default = 3.0,
        )
    nIter = pexConfig.Field(
        dtype = int, 
        doc = 'How many times do we iterate clipping outliers in calculate count statistics?',
        default = 3,
        )

#class QaWriteFits(pexConfig):
class QaConfig(pexConfig.Config):
    seeing = pexConfig.ConfigField(dtype=QaSeeing.QaSeeingConfig, doc="Qa.measSeeing")  
    flatness = pexConfig.ConfigField(dtype=QaFlatness, doc="Qa.flatness")
    doWriteOssImage = pexConfig.Field(
        dtype = bool,
        doc = 'Do we write overscan-subtracted image FITS?',
        default = True,
        )
    doWriteFltImage = pexConfig.Field(
        dtype = bool,
        doc = 'Do we write flatfielded image FITS?',
        default = True,
        )
    doDumpSnapshot = pexConfig.Field(
        dtype=bool,
        doc="Do we dump snapshot files?",
        default=True
        )

class SubaruProcessCcdConfig(ptProcessCcd.ProcessCcdConfig):
    calibrate = pexConfig.ConfigField(dtype=hscCalibrate.HscCalibrateConfig, doc="Calibration")
    qa = pexConfig.ConfigField(dtype=QaConfig, doc="Qa configuration")

class SubaruProcessCcdTask(ptProcessCcd.ProcessCcdTask):
    """Subaru version of ProcessCcdTask, with method to write outputs
    after producing a new multi-frame WCS.
    """
    ConfigClass = SubaruProcessCcdConfig
    
    def write(self, butler, dataId, struct, wcs=None):
        if wcs is None:
            wcs = struct.exposure.getWcs()
            self.log.log(self.log.WARN, "WARNING: No new WCS provided")

        # Apply WCS to sources
        # No longer handling matchSources explicitly - these should all be in calib.sources,
        # or there's a bug in the calibrate task.
        struct.exposure.setWcs(wcs)
        for sources in (struct.sources, struct.calib.sources):
            if sources is None:
                continue
            for s in sources:
                s.updateCoord(wcs)
                
        normalizedMatches = afwTable.packMatches(struct.calib.matches)
        normalizedMatches.table.setMetadata(struct.calib.matchMeta)

        butler.put(struct.exposure, 'calexp', dataId)
        butler.put(struct.sources, 'src', dataId)
        butler.put(normalizedMatches, 'icMatch', dataId)
        butler.put(struct.calib.psf, 'psf', dataId)
        butler.put(struct.calib.apCorr, 'apCorr', dataId)
        butler.put(struct.calib.sources, 'icSrc', dataId)

class SuprimeCamProcessCcdTask(SubaruProcessCcdTask):
    
    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
## FH changed for QA outputs
##        self.makeSubtask("isr", hscSuprimeCam.SuprimeCamIsrTask)
        self.makeSubtask("isr", qaHscSuprimeCamIsrTask.qaSuprimeCamIsrTask)
        self.makeSubtask("calibrate", hscCalibrate.HscCalibrateTask)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                             schema=self.schema, algMetadata=self.algMetadata)

## FH added run() function  for QA output
## derived from & overriding hsc.pipe.tasks.processCcd.ProcessCcdTask(pipeBase.Task).run()
    @pipeBase.timeMethod
    def run(self, sensorRef):
        self.log.log(self.log.INFO, "Processing %s" % (sensorRef.dataId))
        if self.config.doIsr:
            butler = sensorRef.butlerSubset.butler
            calibSet = self.isr.makeCalibDict(butler, sensorRef.dataId)
            rawExposure = sensorRef.get("raw")

            dataId = sensorRef.dataId
##            isrRes = self.isr.run(rawExposure, calibSet)
            isrRes = self.isr.run(rawExposure, calibSet, butler, dataId)            
            self.display("isr", exposure=isrRes.postIsrExposure, pause=True)
            exposure = self.isr.doCcdAssembly([isrRes.postIsrExposure])
            self.display("ccdAssembly", exposure=exposure)
            if self.config.doWriteIsr:
                sensorRef.put(exposure, 'postISRCCD')
        else:
            exposure = None

        if self.config.doCalibrate:
            if exposure is None:
                exposure = sensorRef.get('postISRCCD')
            calib = self.calibrate.run(exposure)
            exposure = calib.exposure
            if self.config.doWriteCalibrate:
                sensorRef.put(exposure, 'calexp')
                sensorRef.put(calib.sources, 'icSrc')
                if calib.psf is not None:
                    sensorRef.put(calib.psf, 'psf')
                if calib.apCorr is not None:
                    sensorRef.put(calib.apCorr, 'apCorr')
                if calib.matches is not None:
                    normalizedMatches = afwTable.packMatches(calib.matches)
                    normalizedMatches.table.setMetadata(calib.matchMeta)
                    sensorRef.put(normalizedMatches, 'icMatch')
        else:
            calib = None

        if self.config.doDetection:
            if exposure is None:
                exposure = sensorRef.get('calexp')
            if calib is None:
                psf = sensorRef.get('psf')
                exposure.setPsf(sensorRef.get('psf'))
            table = afwTable.SourceTable.make(self.schema)
            table.setMetadata(self.algMetadata)
            sources = self.detection.makeSourceCatalog(table, exposure)
        else:
            sources = None

        if self.config.doMeasurement:
            assert(sources)
            assert(exposure)
            if calib is None:
                apCorr = sensorRef.get("apCorr")
            else:
                apCorr = calib.apCorr
            self.measurement.run(exposure, sources, apCorr)
 
        if self.config.doWriteSources:
            sensorRef.put(sources, 'src')

        if self.config.doWriteCalibrate:
            sensorRef.put(exposure, 'calexp')

            
        ##== FH added this part for QA output
        #def measureSeeingQa(exposure, sourceSet, config, debugFlag=False, plotFlag=True, plotbasename=None, io=None, log=None):
        if True:
            QaSeeing.measureSeeingQaTest(exposure, self.config)
        else:
            fwhmRobust, ellRobust, ellPaRobust, sourceSetPsfLike, sourceSetPsfLikeRobust = QaSeeing.measureSeeingQa(exposure, sources, self.config)
            print 'fwhmRobust:', fwhmRobust
            print 'ellRobust:', ellRobust
            print 'ellPaRobust:', ellPaRobust
            print 'sourceSetPsfLike:', sourceSetPsfLike
            print 'sourceSetPsfLikeRobust:', sourceSetPsfLikeRobust
        
        import sys
        sys.exit(0)
        ##==

        return pipeBase.Struct(
            exposure = exposure,
            calib = calib,
            sources = sources,
        )


class HscProcessCcdTask(SubaruProcessCcdTask):

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("isr", IsrTask)
        self.makeSubtask("calibrate", hscDc2.HscDc2CalibrateTask)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                             schema=self.schema, algMetadata=self.algMetadata)
