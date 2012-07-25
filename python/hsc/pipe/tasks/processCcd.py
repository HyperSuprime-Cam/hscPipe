#!/usr/bin/env python

import hsc.pipe.tasks.plotSetup

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
import hsc.pipe.tasks.qa as hscQa

class SubaruProcessCcdConfig(ptProcessCcd.ProcessCcdConfig):
    calibrate = pexConfig.ConfigField(dtype=hscCalibrate.SubaruCalibrateConfig, doc="Calibration")
    isr = pexConfig.ConfigField(dtype=hscIsr.SubaruIsrConfig, doc="Instrument signature removal")
    qa = pexConfig.ConfigField(dtype=hscQa.QaConfig, doc="Quality assessment")
    doWriteUnpackedMatches = pexConfig.Field(
        dtype=bool, default=True,
        doc="Write the denormalized match table as well as the normalized match table"
    )

class SubaruProcessCcdTask(ptProcessCcd.ProcessCcdTask):
    """Subaru version of ProcessCcdTask, with method to write outputs
    after producing a new multi-frame WCS.
    """
    ConfigClass = SubaruProcessCcdConfig
    _DefaultName = "processCcd"

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdTask, self).__init__(*args, **kwargs)
        self.makeSubtask("qa", hscQa.QaTask)

    @pipeBase.timeMethod
    def run(self, sensorRef):
        result = super(SubaruProcessCcdTask, self).run(sensorRef)
        self.qa.run(sensorRef, exposure, sources)
        if self.config.doWriteUnpackedMatches:
            self.writeUnpackedMatches(sensorRef, result.calibrate.matches, result.calibrate.matchMeta)
        return result

    def writeUnpackedMatches(self, dataRef, matches, matchMeta):

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
