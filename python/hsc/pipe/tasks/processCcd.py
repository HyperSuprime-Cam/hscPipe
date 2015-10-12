#!/usr/bin/env python

import os

import hsc.pipe.tasks.plotSetup
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable
import hsc.pipe.base.matches as hscMatches
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

    def setDefaults(self):
        ProcessCcdTask.ConfigClass.setDefaults(self)
        self.doWriteCalibrate = False # We will write it after getting QA metadata


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

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdTask, self).__init__(*args, **kwargs)
        self.makeSubtask("qa")

    def run(self, sensorRef):
        result = ProcessCcdTask.run(self, sensorRef)
        if self.config.qa.useIcsources and self.config.doCalibrate:
            self.qa.run(sensorRef, result.exposure, result.calib.sources)
        elif self.config.doMeasurement:
            self.qa.run(sensorRef, result.exposure, result.sources)

        if self.config.doFinalWrite:
            self.write(sensorRef, result)

        return result

    def write(self, dataRef, results, wcs=None, fluxMag0=None, focusMd=None):
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

        if focusMd is not None:
            results.exposure.getMetadata().combine(focusMd)

        if self.config.doCalibrate:
            if results.calib.matches is not None and results.calib.matchMeta is not None:
                dataRef.put(packMatches(results.calib.matches, results.calib.matchMeta), "icMatch")
                if self.config.doWriteUnpackedMatches:
                    dataRef.put(hscMatches.matchesToCatalog(results.calib.matches, results.calib.matchMeta),
                                "icMatchFull")
            if results.calib.sources is not None:
                dataRef.put(results.calib.sources, 'icSrc')

        if results.backgrounds is not None:
            super(SubaruProcessCcdTask, self).writeBackgrounds(dataRef, results.backgrounds)

        if results.exposure is not None:
            recordCalibInfo(results.exposure, dataRef)
            dataRef.put(results.exposure, 'calexp')
        if results.sources is not None:
            sourceWriteFlags = (0 if self.config.doWriteHeavyFootprintsInSources
                                else afwTable.SOURCE_IO_NO_HEAVY_FOOTPRINTS)
            recordCalibInfo(results.sources, dataRef)
            dataRef.put(results.sources, 'src', flags=sourceWriteFlags)
        if results.matches is not None:
            dataRef.put(packMatches(results.matches, results.matchMeta), "srcMatch")
            if self.config.doWriteUnpackedMatches:
                dataRef.put(hscMatches.matchesToCatalog(results.matches, results.matchMeta), "srcMatchFull")

def recordCalibInfo(obj, dataRef):
    """Local function to add some header information about calibration inputs to the detrend outputs

    This could go in butler.put(), but that fails for detrend creation (all calibVersions are set
    to the output calibVersion).  So we'll do this in processCcd, and then use a duplicated separate
    one in detrend.py
    """

    try:
        md = obj.getMetadata()
    except AttributeError:
        return
    if md is None:
        return

    for calibName in 'bias', 'dark', 'flat', 'fringe':

        calibPath, calibVersion, calibDate = "not_available", "not_available", "not_available"
        try:
            calibPath      = dataRef.get(calibName+"_filename")
            additionalData = dataRef.getButler().mapper.map(calibName, dataRef.dataId).getAdditionalData()
            calibVersion   = additionalData.get('calibVersion')
            calibDate      = additionalData.get('calibDate')
        except:
            pass
        md.add(calibName.upper()+"_VERSION", calibVersion)
        md.add(calibName.upper()+"_DATE", calibDate)
        md.add(calibName.upper()+"_PATH", calibPath)
