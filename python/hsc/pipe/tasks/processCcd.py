#!/usr/bin/env python

import os

import hsc.pipe.tasks.plotSetup
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
from hsc.pipe.base import SubaruArgumentParser
from lsst.pipe.tasks.processCcd import ProcessCcdTask
from .qa import QaTask


class SubaruProcessCcdConfig(ProcessCcdTask.ConfigClass):
    delayWrite = pexConfig.Field(
        dtype=bool, default=False,
        doc="Delay writing outputs (e.g., for pbasf processing)?"
        )
    doWriteUnpackedMatches = pexConfig.Field(
        dtype=bool, default=True,
        doc=("Write the denormalized match table as well as the normalized match table; "
             "ignored if doWriteCalibrate=False")
    )
    qa = pexConfig.ConfigurableField(target = QaTask, doc = "QA analysis")


def applyOverrides(root):
    root.doWriteCalibrate = False # We will write it after getting QA metadata


class SubaruProcessCcdTask(ProcessCcdTask):
    """Subaru version of ProcessCcdTask, with method to write outputs
    after producing a new multi-frame WCS, the ability to write denormalized
    matches, and --rerun support.
    """
    ConfigClass = SubaruProcessCcdConfig
    overrides = (applyOverrides,)

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdTask, self).__init__(*args, **kwargs)
        self.makeSubtask("qa")

    def run(self, sensorRef):
        result = ProcessCcdTask.run(self, sensorRef)
        self.qa.run(sensorRef, result.exposure, result.sources)
        sensorRef.put(result.exposure, self.dataPrefix + 'calexp')
        if not self.config.delayWrite:
            self.write(sensorRef, result, wcs=result.exposure.getWcs())
        if self.config.doWriteUnpackedMatches:
            sensorRef.put(self.unpackMatches(result.calib.matches, result.calib.matchMeta), "icMatchList")
        if self.config.doWriteSourceMatches and self.config.doWriteUnpackedMatches:
            sensorRef.put(self.unpackMatches(result.matches, result.matchMeta), "srcMatchList")

        return result

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser with --rerun support.
        """
        return SubaruArgumentParser(name=cls._DefaultName)

    def unpackMatches(self, matches, matchMeta):
        """Denormalise matches into "unpacked matches" """

        refSchema = matches[0].first.getSchema()
        srcSchema = matches[0].second.getSchema()

        mergedSchema = afwTable.Schema()
        def merge(schema, key, merged, name):
            field = schema.find(key).field
            typeStr = field.getTypeString()
            fieldDoc = field.getDoc()
            fieldUnits = field.getUnits()
            if typeStr in ("ArrayF", "ArrayD", "CovF", "CovD"):
                fieldSize = field.getSize()
            else:
                fieldSize = None
            mergedSchema.addField(name + '.' + key, type=typeStr, doc=fieldDoc, units=fieldUnits,
                                  size=fieldSize)

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

        for match in matches:
            record = mergedCatalog.addNew()
            for key in refKeys:
                keyIn = key[0]
                keyOut = key[1]
                record.set(keyOut, match.first.get(keyIn))
            for key in srcKeys:
                keyIn = key[0]
                keyOut = key[1]
                record.set(keyOut, match.second.get(keyIn))

        # obtaining reference catalog name
        catalogName = os.path.basename(os.getenv("ASTROMETRY_NET_DATA_DIR").rstrip('/'))
        matchMeta.add('REFCAT', catalogName)
        mergedCatalog.getTable().setMetadata(matchMeta)

        return mergedCatalog

    def write(self, dataRef, struct, wcs=None):
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
        dataRef.put(self.unpackMatches(struct.calib.matches, struct.calib.matchMeta), "icMatchList")
        dataRef.put(struct.exposure, 'calexp')
        dataRef.put(struct.sources, 'src')
        dataRef.put(normalizedMatches, "icMatch")
        dataRef.put(struct.calib.psf, 'psf')
        dataRef.put(struct.calib.apCorr, 'apCorr')
        dataRef.put(struct.calib.sources, 'icSrc')
