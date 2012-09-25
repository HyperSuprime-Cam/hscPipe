#!/usr/bin/env python

import os

import hsc.pipe.tasks.plotSetup
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import hsc.pipe.base.matches as hscMatches
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
        if self.config.qa.useIcsources:
            self.qa.run(sensorRef, result.exposure, result.calib.sources)        
        else:
            self.qa.run(sensorRef, result.exposure, result.sources)
        sensorRef.put(result.exposure, self.dataPrefix + 'calexp')
        if not self.config.delayWrite:
            self.write(sensorRef, result, wcs=result.exposure.getWcs())
        if self.config.doWriteUnpackedMatches:
            sensorRef.put(hscMatches.matchesToCatalog(result.calib.matches, result.calib.matchMeta), "icMatchList")
        if self.config.doWriteSourceMatches and self.config.doWriteUnpackedMatches and (not self.config.qa.useIcsources):
            sensorRef.put(hscMatches.matchesToCatalog(result.matches, result.matchMeta), "srcMatchList")

        return result

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser with --rerun support.
        """
        return SubaruArgumentParser(name=cls._DefaultName)

    def write(self, dataRef, struct, wcs=None, fluxMag0=None):
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

        if fluxMag0 is not None:
            struct.exposure.getCalib().setFluxMag0(fluxMag0)

        if self.config.doCalibrate and self.config.doWriteCalibrate:
            if struct.calib.psf is not None:
                dataRef.put(struct.calib.psf, 'psf')
            if self.config.doWriteCalibrateMatches and struct.calib.matches is not None and struct.calib.matchMeta is not None:
                normalizedMatches = afwTable.packMatches(struct.calib.matches)
                normalizedMatches.table.setMetadata(struct.calib.matchMeta)
                dataRef.put(normalizedMatches, "icMatch")
            if self.config.doWriteUnpackedMatches and  struct.calib.matches is not None and struct.calib.matchMeta is not None:
                dataRef.put(hscMatches.matchesToCatalog(struct.calib.matches, struct.calib.matchMeta),
                            "icMatchList")
            if self.config.calibrate.doComputeApCorr and struct.calib.apCorr is not None:
                dataRef.put(struct.calib.apCorr, 'apCorr')
            if struct.calib.sources is not None:
                dataRef.put(struct.calib.sources, 'icSrc')

        if struct.exposure is not None:
            dataRef.put(struct.exposure, 'calexp')
        if self.config.doWriteSources and not self.config.qa.useIcsources and struct.sources is not None:
            dataRef.put(struct.sources, 'src')
