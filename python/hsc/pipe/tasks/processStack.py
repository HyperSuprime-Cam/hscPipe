#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import os
import numpy

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import hsc.pipe.base as hscPipeBase
import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
from lsst.pipe.tasks.processImage import ProcessImageTask

class ProcessStackConfig(ProcessImageTask.ConfigClass):
    """Config for ProcessStack"""
    coaddName = pexConfig.Field(
        doc = "coadd name: typically one of deep or goodSeeing",
        dtype = str,
        default = "deep",
    )
    doScaleVariance = pexConfig.Field(dtype=bool, default=True,
                                      doc = "Scale variance plane using empirical noise")
    #delayWrite = pexConfig.Field(
    #            dtype=bool, default=False,
    #                    doc="Delay writing outputs (e.g., for pbasf processing)?"
    #                    )
    doWriteUnpackedMatches = pexConfig.Field(
        dtype=bool, default=True,
        doc=("Write the denormalized match table as well as the normalized match table; "
             "ignored if doWriteCalibrate=False")
        )

class ProcessStackTask(ProcessImageTask):
    """Process a Stack image
    
    """
    ConfigClass = ProcessStackConfig
    _DefaultName = "processStack"

    def __init__(self, **kwargs):
        ProcessImageTask.__init__(self, **kwargs)
        self.dataPrefix = "stack_"

    @pipeBase.timeMethod
    def scaleVariance(self, exposure):
        ctrl = afwMath.StatisticsControl()
        ctrl.setAndMask(~0x0)
        var    = exposure.getMaskedImage().getVariance()
        mask   = exposure.getMaskedImage().getMask()
        dstats = afwMath.makeStatistics(exposure.getMaskedImage(), afwMath.VARIANCECLIP, ctrl).getValue(afwMath.VARIANCECLIP)
        vstats = afwMath.makeStatistics(var, mask, afwMath.MEANCLIP, ctrl).getValue(afwMath.MEANCLIP)
        vrat   = dstats / vstats
        self.log.info("Renormalising variance by %f" % (vrat))
        var   *= vrat

    @pipeBase.timeMethod
    def fixEdgeNaNs(self, exposure):
        """Replace NaN pixels that are also marked as EDGE with 0.

        It might be more efficient to do this with C++ code, but NumPy makes it
        really easy to write.
        """
        self.log.log(self.log.INFO, "Setting NaN EDGE pixels to zero.")
        image = exposure.getMaskedImage().getImage().getArray()
        mask = exposure.getMaskedImage().getMask().getArray()
        edgeBitMask = afwImage.MaskU.getPlaneBitMask("EDGE")
        toReplace = numpy.logical_and(
            numpy.isnan(image),
            mask & edgeBitMask
            )
        image[toReplace] = 0.0

    def makeIdFactory(self, dataRef):
        expBits = dataRef.get("stackExposureId_bits")
        expId = long(dataRef.get("stackExposureId"))
        return afwTable.IdFactory.makeSource(expId, 64 - expBits)

    @pipeBase.timeMethod
    def run(self, dataRef):
        """Process a coadd image
        
        @param dataRef: butler data reference corresponding to coadd patch
        @return pipe_base Struct containing these fields:
        - exposure: calibrated exposure (calexp): as computed if config.doCalibrate,
            else as upersisted and updated if config.doDetection, else None
        - calib: object returned by calibration process if config.doCalibrate, else None
        - apCorr: aperture correction: as computed config.doCalibrate, else as unpersisted
            if config.doMeasure, else None
        - sources: detected source if config.doDetection, else None
        """
        self.log.log(self.log.INFO, "Processing %s" % (dataRef.dataId))

        # initialize outputs
        calExposure = None
        calib = None
        apCorr = None
        psf = None

        if self.config.doCalibrate:
            coadd = dataRef.get("stack")
            if dataRef.datasetExists(self.dataPrefix + "initPsf"):
                initPsf = dataRef.get(self.dataPrefix + "initPsf")
                coadd.setPsf(initPsf)
            else:
                self.log.warn("Could not load initial PSF; dataset does not exist")
            if self.config.doScaleVariance:
                self.scaleVariance(coadd)
            self.fixEdgeNaNs(coadd)
        else:
            coadd = None

        # delegate most of the work to ProcessImageTask
        result = self.process(dataRef, coadd)
        result.coadd = coadd

        if self.config.doWriteUnpackedMatches:
            dataRef.put(self.unpackMatches(result.calib.matches, result.calib.matchMeta), self.dataPrefix + "icMatchFull")
        if self.config.doWriteSourceMatches and self.config.doWriteUnpackedMatches:
            dataRef.put(self.unpackMatches(result.matches, result.matchMeta), self.dataPrefix + "srcMatchFull")

        return result


    # The below unpackMatches function is a wholesale copy from processCcd
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

    @classmethod
    def _makeArgumentParser(cls):
        return StackArgumentParser(name=cls._DefaultName, datasetType="stack", defaultToUserName=False)


class StackArgumentParser(hscPipeBase.SubaruArgumentParser):
    """A version of lsst.pipe.base.ArgumentParser specialized for stacks.
    
    Required because butler.subset does not support patch and tract
    """

    def _makeDataRefList(self, namespace):
        """Make namespace.dataRefList from namespace.dataIdList
        """
        validKeys = namespace.butler.getKeys(datasetType="stack", level=self._dataRefLevel)

        namespace.dataRefList = []
        for dataId in namespace.dataIdList:
            # stack, patch, and filter are required
            for key in validKeys:
                if key not in dataId:
                    self.error("--id must include " + key)
            dataRef = namespace.butler.dataRef(
                datasetType = "stack",
                dataId = dataId,
            )
            namespace.dataRefList.append(dataRef)

