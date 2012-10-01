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
    doFinalWrite = pexConfig.Field(
        dtype=bool, default=True,
        doc="Write all outputs at end?  If so, you also likely want doWriteCalibrate=False."
        )
    doWriteUnpackedMatches = pexConfig.Field(
        dtype=bool, default=True,
        doc=("Write the denormalized match table as well as the normalized match table (in the final write)?")
    )
    qa = pexConfig.ConfigurableField(target = QaTask, doc = "QA analysis")


def applyOverrides(root):
    root.doWriteCalibrate = False # We will write it after getting QA metadata


def packMatches(matches, matchMeta):
    normalizedMatches = afwTable.packMatches(matches)
    normalizedMatches.table.setMetadata(matchMeta)
    return normalizedMatches

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

        if self.config.doFinalWrite:
            self.write(sensorRef, result, wcs=result.exposure.getWcs())

        return result

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser with --rerun support.
        """
        return SubaruArgumentParser(name=cls._DefaultName)

    def write(self, dataRef, results, wcs=None, fluxMag0=None):
        if wcs is None:
            self.log.warn("WARNING: No new WCS provided")
        else:
            # Apply WCS to sources
            # No longer handling matchSources explicitly - these should all be in calib.sources,
            # or there's a bug in the calibrate task.
            if results.exposure is not None:
                results.exposure.setWcs(wcs)
            for sources in (results.sources, results.calib.sources):
                if sources is None:
                    continue
                for s in sources:
                    s.updateCoord(wcs)

        if fluxMag0 is not None:
            results.exposure.getCalib().setFluxMag0(fluxMag0)

        if self.config.doCalibrate:
            if results.calib.psf is not None:
                dataRef.put(results.calib.psf, 'psf')
            if results.calib.matches is not None and results.calib.matchMeta is not None:
                dataRef.put(packMatches(results.calib.matches, results.calib.matchMeta), "icMatch")
                if self.config.doWriteUnpackedMatches:
                    dataRef.put(hscMatches.matchesToCatalog(results.calib.matches, results.calib.matchMeta),
                                "icMatchFull")
            if self.config.calibrate.doComputeApCorr and results.calib.apCorr is not None:
                dataRef.put(results.calib.apCorr, 'apCorr')
            if results.calib.sources is not None:
                dataRef.put(results.calib.sources, 'icSrc')

        if results.exposure is not None:
            dataRef.put(results.exposure, 'calexp')
        if results.sources is not None:
            dataRef.put(results.sources, 'src')
        if results.matches is not None:
            dataRef.put(packMatches(results.matches, results.matchMeta), "srcMatch")
            if self.config.doWriteUnpackedMatches:
                dataRef.put(hscMatches.matchesToCatalog(results.matches, results.matchMeta), "srcMatchFull")
