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
import hsc.pipe.tasks.qaSuprimeCamIsrTask as qaSuprimeCamIsrTask
import hsc.pipe.tasks.calibrate as hscCalibrate
import hsc.pipe.tasks.hscDc2 as hscDc2

##== FH added for QA output
import hsc.onsite.qa.measSeeingQa as qaSeeing
import hsc.pipe.tasks.qaCalibrate as qaHscCalibrate

#== FH changed for QA output

#class QaWriteFits(pexConfig):
class QaConfig(pexConfig.Config):
    seeing = pexConfig.ConfigField(dtype=qaSeeing.QaSeeingConfig, doc="Qa.measSeeing")  

class SubaruProcessCcdConfig(ptProcessCcd.ProcessCcdConfig):
    calibrate = pexConfig.ConfigField(dtype=qaHscCalibrate.HscCalibrateConfig, doc="Calibration")
##    calibrate = pexConfig.ConfigField(dtype=hscCalibrate.HscCalibrateConfig, doc="Calibration")    
    isr = pexConfig.ConfigField(dtype=qaSuprimeCamIsrTask.QaIsrTaskConfig, doc="Config for Isr with Qa tasks")
    qa = pexConfig.ConfigField(dtype=QaConfig, doc="Config for Qa (outside isr) tasks")

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
##        self.makeSubtask("calibrate", hscCalibrate.HscCalibrateTask)
        self.makeSubtask("isr", qaSuprimeCamIsrTask.QaSuprimeCamIsrTask)#, config=SubaruProcessCcdConfig())
        self.makeSubtask("calibrate", qaHscCalibrate.HscCalibrateTask) #, config=self.config)        
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
                    print '*** len(normalizedMatches):',len(normalizedMatches)
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

        ##== FH added this part for QA output by using the reduced exposure and final sources
        ## debug: QaSeeing.measureSeeingQaTest(exposure, self.config)
        fwhmRobust, ellRobust, ellPaRobust, catalogPsfLike, catalogPsfLikeRobust = qaSeeing.measureSeeingQa(exposure, sources, self.config, debugFlag=False, plotFlag=True, plotbasedir=None, butler=butler, log=self.log)

        self.log.log(self.log.INFO, "QA seeing: fwhm: %f (pix)" % fwhmRobust)
        self.log.log(self.log.INFO, "QA seeing: ell (based on 2nd moments): %f" % ellRobust)
        self.log.log(self.log.INFO, "QA seeing: ellPa (in CCDCoords based on 2nd moments): %f (deg)" % ellRobust)
        self.log.log(self.log.INFO, "QA seeing: final Nsources for seeing: %d" % len(catalogPsfLikeRobust))        

        # this part should be done by calculating merit functions somewhere else in a polite manner.
        metadata = exposure.getMetadata()
        metadata.set('FLAG_AUTO', 0)
        metadata.set('FLAG_USR', 0)
        metadata.set('FLAG_TAG', 1)

 
        if self.config.doWriteSources:
            sensorRef.put(sources, 'src')

        if self.config.doWriteCalibrate:
            sensorRef.put(exposure, 'calexp')

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
