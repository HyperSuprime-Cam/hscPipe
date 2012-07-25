#!/usr/bin/env python

import os, os.path
from lsst.pex.config import Field
from lsst.ip.isr.isrTask import IsrTaskConfig, IsrTask
from lsst.ip.isr.isr import Isr
from lsst.pipe.base import Struct

import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import hsc.fitsthumb as fitsthumb

import numpy

class SubaruIsrConfig(IsrTaskConfig):
    doSaturation = Field(doc="Mask saturated pixels?", dtype=bool, default=True)
    doOverscan = Field(doc="Do overscan subtraction?", dtype=bool, default=True)
    doBias = Field(doc="Do bias subtraction?", dtype=bool, default=False)
    doVariance = Field(doc="Calculate variance?", dtype=bool, default=True)
    doDark = Field(doc="Do dark subtraction?", dtype=bool, default=False)
    doFlat = Field(doc="Do flat-fielding?", dtype=bool, default=True)
    doWriteOss = Field(doc="Write OverScan-Subtracted image?", dtype=bool, default=False)
    doThumbnailOss = Field(doc="Write OverScan-Subtracted thumbnail?", dtype=bool, default=True)
    doWriteFlattened = Field(doc="Write flattened image?", dtype=bool, default=False)
    doThumbnailFlattened = Field(doc="Write flattened thumbnail?", dtype=bool, default=True)
    meshX = Field(dtype=int, doc='Mesh size in X (pix) to calculate count statistics', default=256)
    meshY = Field(dtype=int, doc='Mesh size in Y (pix) to calculate count statistics', default=256)
    doClip = Field(dtype=bool, doc='Do we clip outliers in calculate count statistics?', default=True)
    clipSigma = Field(dtype=float, doc='Clipping threshold (sigma)', default=3.0)
    nIter = Field(dtype=int, doc='Clipping iterations', default=3)

class SubaruIsr(Isr):
    def overscanCorrection(self, maskedImage, overscanData, *args, **kwargs):
        stats = afwMath.makeStatistics(overscanData, afwMath.MEDIAN|afwMath.STDEVCLIP)
        osLevel = stats.getValue(afwMath.MEDIAN)
        osSigma = stats.getValue(afwMath.STDEVCLIP)

        super(SubaruIsr, self).overscanCorrection(maskedImage, overscanData, *args, **kwargs)
        return 


class SubaruIsrTask(IsrTask):
    ConfigClass = SubaruIsrConfig
    def __init__(self, *args, **kwargs):
        super(SubaruIsrTask, self).__init__(*args, **kwargs)
        self.isr = SubaruIsr()

    def run(self, dataRef, exposure, calibSet):
        exposure = self.doConversionForIsr(exposure, calibSet)

        self.measureOverscan(exposure)

        if self.config.doSaturation:
            exposure = self.doSaturationDetection(exposure, calibSet)
        if self.config.doOverscan:
            exposure = self.doOverscanCorrection(exposure, calibSet)

        if self.config.doVariance:
            # Ideally, this should be done after bias subtraction, but CCD assembly demands a variance plane
            exposure = self.doVariance(exposure, calibSet)

        exposure = self.doCcdAssembly([exposure])

        if self.config.doWriteOss:
            dataRef.put("ossImage", exposure)
        if self.config.doThumbnailOss:
            self.writeThumbnail(dataRef, "ossThumb", exposure)

        if self.config.doBias:
            exposure = self.doBiasSubtraction(exposure, calibSet)
        if self.config.doDark:
            exposure = self.doDarkCorrection(exposure, calibSet)
        if self.config.doFlat:
            exposure = self.doFlatCorrection(exposure, calibSet)

        if self.config.doWriteFlattened:
            dataRef.put("flattenedImage", exposure)
        if self.config.doThumbnailFlattened:
            self.writeThumbnail(dataRef, "flattenedThumb", exposure)

        self.measureBackground(exposure)

        return Struct(postIsrExposure=exposure)

    def makeCalibDict(self, butler, dataId):
        ret = {}
        required = {"doBias": "bias",
                    "doDark": "dark",
                    "doFlat": "flat",
                    }
        for method in required.keys():
            if getattr(self.config, method):
                calib = required[method]
                ret[calib] = butler.get(calib, dataId)
        return ret


    def writeThumbnail(self, dataRef, dataset, exposure, format='png', width=500, height=0):
        """Write out exposure to a snapshot file named outfile in the given image format and size.
        """
        filename = dataRef.get(dataset + "_filename")[0]
        directory = os.path.dirname(filename)
        if not os.path.exists(directory):
            os.makedirs(directory)
        image = exposure.getMaskedImage().getImage()
        fitsthumb.createFitsThumb(image, filename, format, width, height, True)

    def measureOverscan(self, exposure):
        clipSigma = 3.0
        nIter = 3
        levelStat = afwMath.MEDIAN
        sigmaStat = afwMath.STDEVCLIP
        
        metadata = exposure.getMetadata()
        sctrl = afwMath.StatisticsControl(clipSigma, nIter)
        for amp in self._getAmplifiers(exposure):
            expImage = exposure.getMaskedImage().getImage()
            overscan = expImage.Factory(expImage, amp.getDiskBiasSec())
            stats = afwMath.makeStatistics(overscan, levelStat | sigmaStat, sctrl)
            ampNum = amp.getId().getSerial()
            metadata.set("OSLEVEL%d" % ampNum, stats.getValue(levelStat))
            metadata.set("OSSIGMA%d" % ampNum, stats.getValue(sigmaStat))


    def measureBackground(self, exposure):
        statsControl = afwMath.StatisticsControl(self.config.clipSigma, self.config.nIter)
        maskedImage = exposure.getMaskedImage()
        stats = afwMath.makeStatistics(maskedImage, afwMath.MEDIAN | afwMath.STDEVCLIP, statsControl)
        skyLevel = stats.getValue(afwMath.MEDIAN)
        skySigma = stats.getValue(afwMath.STDEVCLIP)
        self.log.info("Flattened sky level: %f +/- %f" % (skyLevel, skySigma))
        metadata = exposure.getMetadata()
        metadata.set('SKYLEVEL', skyLevel)
        metadata.set('SKYSIGMA', skySigma)

        # calcluating flatlevel over the subgrids 
        stat = afwMath.MEANCLIP if self.config.doClip else afwMath.MEAN
        meshXHalf = int(self.config.meshX/2.)
        meshYHalf = int(self.config.meshY/2.)
        nX = int((exposure.getWidth() + meshXHalf) / self.config.meshX)
        nY = int((exposure.getHeight() + meshYHalf) / self.config.meshY)
        skyLevels = numpy.zeros((nX,nY))

        for j in range(nY):
            yc = meshYHalf + j * self.config.meshY
            for i in range(nX):
                xc = meshXHalf + i * self.config.meshX

                xLLC = xc - meshXHalf
                yLLC = yc - meshYHalf
                xURC = xc + meshXHalf - 1
                yURC = yc + meshYHalf - 1

                bbox = afwGeom.Box2I(afwGeom.Point2I(xLLC, yLLC), afwGeom.Point2I(xURC, yURC))
                miMesh = maskedImage.Factory(exposure.getMaskedImage(), bbox, afwImage.LOCAL)

                skyLevels[i,j] = afwMath.makeStatistics(miMesh, stat, statsControl).getValue()

        skyMedian = numpy.median(skyLevels)
        flatness =  (skyLevels - skyMedian) / skyMedian
        flatness_rms = numpy.std(flatness)
        flatness_min = flatness.min()
        flatness_max = flatness.max() 
        flatness_pp = flatness_max - flatness_min

        self.log.info("Measuring sky levels in %dx%d grids: %f" % (nX, nY, skyMedian))
        self.log.info("Sky flatness in %dx%d grids - pp: %f rms: %f" % (nX, nY, flatness_pp, flatness_rms))

        metadata.set('FLATNESS_PP', flatness_pp)
        metadata.set('FLATNESS_RMS', flatness_rms)
        metadata.set('FLATNESS_NGRIDS', '%dx%d' % (nX, nY))
        metadata.set('FLATNESS_MESHX', self.config.meshX)
        metadata.set('FLATNESS_MESHY', self.config.meshY)

