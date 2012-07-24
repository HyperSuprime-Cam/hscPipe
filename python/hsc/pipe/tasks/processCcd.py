#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDet
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.processCcd as ptProcessCcd
import hsc.pipe.tasks.suprimecam as hscSuprimeCam
import hsc.pipe.tasks.calibrate as hscCalibrate
import hsc.pipe.tasks.isr as hscIsr
import hsc.pipe.tasks.hscDc2 as hscDc2
import hsc.pipe.tasks.qa as hscQa


class SubaruProcessCcdConfig(ptProcessCcd.ProcessCcdConfig):
    calibrate = pexConfig.ConfigField(dtype=hscCalibrate.SubaruCalibrateConfig, doc="Calibration")
    isr = pexConfig.ConfigField(dtype=hscIsr.SubaruIsrConfig, doc="Instrument signature removal")
    qa = pexConfig.ConfigField(dtype=hscQa.QaConfig, doc="Quality assessment")


class SubaruProcessCcdTask(ptProcessCcd.ProcessCcdTask):
    """Subaru version of ProcessCcdTask, with method to write outputs
    after producing a new multi-frame WCS.
    """
    ConfigClass = SubaruProcessCcdConfig
    _DefaultName = "processCcd"

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdTask, self).__init__(*args, **kwargs)
        self.makeSubtask("isr", hscIsr.SubaruIsrTask)
        self.makeSubtask("calibrate", hscCalibrate.SubaruCalibrateTask)
        self.makeSubtask("qa", hscQa.QaTask)

    # The 'run' method is copied wholesale from lsst.pipe.tasks.processCcd.ProcessCcdTask.run, with minor
    # modifications to change when the CCD assembly is performed.
    @pipeBase.timeMethod
    def run(self, sensorRef):
        self.log.log(self.log.INFO, "Processing %s" % (sensorRef.dataId))
        if self.config.doIsr:
            butler = sensorRef.butlerSubset.butler
            calibSet = self.isr.makeCalibDict(butler, sensorRef.dataId)
            exposure = sensorRef.get("raw")
            isrRes = self.isr.run(sensorRef, exposure, calibSet)
            exposure = isrRes.postIsrExposure
            self.display("isr", exposure=exposure, pause=True)
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
            detRet = self.detection.makeSourceCatalog(table, exposure)
            sources = detRet.sources
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

        return pipeBase.Struct(
            exposure = exposure,
            calib = calib,
            sources = sources,
        )

    def writeCalib(self, sensorRef, calib):
        """Write calibration products
       
        @param sensorRef   Data reference for sensor
        @param calib       Results of calibration
        """
        sensorRef.put(calib.sources, 'icSrc')
        if calib.psf is not None:
            sensorRef.put(calib.psf, 'psf')
        if calib.apCorr is not None:
            sensorRef.put(calib.apCorr, 'apCorr')
        if calib.matches is not None:
            self.writeMatches(sensorRef, calib.matches, calib.matchMeta)

    def writeMatches(self, dataRef, matches, matchMeta):
        # First write normalised matches
        normalizedMatches = afwTable.packMatches(calib.matches)
        normalizedMatches.table.setMetadata(calib.matchMeta)
        dataRef.put(normalizedMatches, 'icMatch')

        # Now write unpacked matches
        refSchema = matchlist[0].first.getSchema()
        srcSchema = matchlist[0].second.getSchema()

        mergedSchema = afwTable.Schema()
        def merge(schema, key, merged, name):
            field = schema.find(key).field
            typeStr = field.getTypeString()
            fieldDoc = field.getDoc()
            fieldUnits = field.getUnits()
            mergedSchema.addField(name + '.' + key, type=typeStr, doc=fieldDoc, units=fieldUnits)

        for keyName in refSchema.getNames():
            merge(refSchema, keyName, mergedSchema, "ref")

        for keyName in srcSchema.getNames():
            merge(srcSchema, keyName, mergedSchema, "src")

        mergedCatalog = afwTable.BaseCatalog(mergedSchema)

        refKeys = []
        for keyName in refSchema.getNames():
            refKeys.append((refSchema.find(keyName).key, mergedSchema.find('ref.' + keyName).key))
        srcKeys = []
        for keyName in srcSchema.getNames():
            srcKeys.append((srcSchema.find(keyName).key, mergedSchema.find('src.' + keyName).key))

        # obtaining reference catalog name
        catalogName = os.path.basename(os.getenv("ASTROMETRY_NET_DATA_DIR").rstrip('/'))
        matchMeta.add('REFCAT', catalogName)
        mergedCatalog.getTable().setMetadata(matchMeta)

        dataRef.put(mergedCatalog, "matchedList")


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

        self.writeMatches(struct.calib.matches, struct.calib.matchMeta)

        butler.put(struct.exposure, 'calexp', dataId)
        butler.put(struct.sources, 'src', dataId)
        butler.put(normalizedMatches, 'icMatch', dataId)
        butler.put(struct.calib.psf, 'psf', dataId)
        butler.put(struct.calib.apCorr, 'apCorr', dataId)
        butler.put(struct.calib.sources, 'icSrc', dataId)


class SuprimeCamProcessCcdConfig(SubaruProcessCcdConfig):
    isr = pexConfig.ConfigField(dtype=hscSuprimeCam.SuprimeCamIsrTask.ConfigClass,
                                doc="Instrument signature removal")

class SuprimeCamProcessCcdTask(SubaruProcessCcdTask):
    ConfigClass = SuprimeCamProcessCcdConfig

    def __init__(self, **kwargs):
        pipeBase.Task.__init__(self, **kwargs)
        self.makeSubtask("isr", hscSuprimeCam.SuprimeCamIsrTask)
        self.makeSubtask("calibrate", hscCalibrate.SubaruCalibrateTask)
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                             schema=self.schema, algMetadata=self.algMetadata)


class HscProcessCcdTask(SubaruProcessCcdTask):
   def __init__(self, **kwargs):
       pipeBase.Task.__init__(self, **kwargs)
       self.makeSubtask("isr", hscIsr.SubaruIsrTask)
       self.makeSubtask("calibrate", hscCalibrate.SubaruCalibrateTask)
       self.schema = afwTable.SourceTable.makeMinimalSchema()
       self.algMetadata = dafBase.PropertyList()
       if self.config.doDetection:
           self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
       if self.config.doMeasurement:
           self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                            schema=self.schema, algMetadata=self.algMetadata)

