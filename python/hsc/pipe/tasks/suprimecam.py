#!/usr/bin/env python

import numpy as np
import lsst.pex.config as pexConfig
import lsst.afw.cameraGeom as afwCG
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.ip.isr as ipIsr
import lsst.pipe.tasks.processCcd as ptProcessCcd
import hsc.pipe.tasks.calibrate as hscCalibrate
import hsc.pipe.tasks.isr as hscIsr
import lsst.meas.algorithms.crosstalk as maCrosstalk

class HamamatsuIsrTaskConfig(hscIsr.HscIsrTask.ConfigClass):
    doCrosstalk = pexConfig.Field(
        dtype = bool,
        doc = "Correct for crosstalk",
        default = True,
    )
    crosstalkCoeffs = pexConfig.ConfigField(
        dtype = maCrosstalk.CrosstalkCoeffsConfig,
        doc = "Crosstalk coefficients",
    )

    crosstalkMaskPlane = pexConfig.Field(
        dtype = str,
        doc = "Name of Mask plane for crosstalk corrected pixels",
        default = "CROSSTALK",
    )

    minPixelToMask = pexConfig.Field(
        dtype = float,
        doc = "Minimum pixel value (in electrons) to cause crosstalkMaskPlane bit to be set",
        default = 45000,
        )

    doLinearize = pexConfig.Field(
        dtype = bool,
        doc = "Correct for nonlinearity of the detector's response (ignored if coefficients are 0.0)",
        default = True,
    )
    linearizationThreshold = pexConfig.Field(
        dtype = float,
        doc = "Minimum pixel value (in electrons) to apply linearity corrections",
        default = 0.0,
        )
    linearizationCoefficient = pexConfig.Field(
        dtype = float,
        doc = "Linearity correction coefficient",
        default = 0.0,
        )

class SuprimeCamIsrTaskConfig(HamamatsuIsrTaskConfig):
    pass

class SuprimeCamIsrTask(hscIsr.HscIsrTask):
    ConfigClass = SuprimeCamIsrTaskConfig

    def run(self, *args, **kwargs):
        result = super(SuprimeCamIsrTask, self).run(*args, **kwargs)
        exposure = result.postIsrExposure

        if self.config.doCrosstalk:
            self.crosstalk(exposure)

        if self.config.doLinearize:
            self.linearize(exposure)

        self.guider(exposure)
        return result

    def crosstalk(self, exposure):
        coeffs = self.config.crosstalkCoeffs.getCoeffs()

        if np.any(coeffs):
            self.log.log(self.log.INFO, "Applying crosstalk corrections to CCD %s" %
                         (exposure.getDetector().getId()))

            maCrosstalk.subtractXTalk(exposure.getMaskedImage(), coeffs,
                                      self.config.minPixelToMask, self.config.crosstalkMaskPlane)

    def guider(self, exposure):
        """Mask defects and trim guider shadow

        @param exposure Exposure to process
        @return Defect list
        """
        assert exposure, "No exposure provided"

        ccd = afwCG.cast_Ccd(exposure.getDetector()) # This is Suprime-Cam so we know the Detector is a Ccd
        ccdNum = ccd.getId().getSerial()
        if ccdNum not in [0, 1, 2, 6, 7]:
            # No need to mask
            return

        md = exposure.getMetadata()
        if not md.exists("S_AG-X"):
            self.log.log(self.log.WARN, "No autoguider position in exposure metadata.")
            return

        xGuider = md.get("S_AG-X")
        if ccdNum in [1, 2, 7]:
            maskLimit = int(60.0 * xGuider - 2300.0) # From SDFRED
        elif ccdNum in [0, 6]:
            maskLimit = int(60.0 * xGuider - 2000.0) # From SDFRED

        mi = exposure.getMaskedImage()
        height = mi.getHeight()
        if height < maskLimit:
            # Nothing to mask!
            return

        if False:
            # XXX This mask plane isn't respected by background subtraction or source detection or measurement
            self.log.log(self.log.INFO, "Masking autoguider shadow at y > %d" % maskLimit)
            mask = mi.getMask()
            bbox = afwGeom.Box2I(afwGeom.Point2I(0, maskLimit - 1),
                                 afwGeom.Point2I(mask.getWidth() - 1, height - 1))
            badMask = mask.Factory(mask, bbox, afwImage.LOCAL)

            mask.addMaskPlane("GUIDER")
            badBitmask = mask.getPlaneBitMask("GUIDER")

            badMask |= badBitmask
        else:
            # XXX Temporary solution until a mask plane is respected by downstream processes
            self.log.log(self.log.INFO, "Removing pixels affected by autoguider shadow at y > %d" % maskLimit)
            bbox = afwGeom.Box2I(afwGeom.Point2I(0, 0), afwGeom.Extent2I(mi.getWidth(), maskLimit))
            good = mi.Factory(mi, bbox, afwImage.LOCAL)
            exposure.setMaskedImage(good)

    def linearize(self, exposure):
        """Correct for non-linearity

        @param exposure Exposure to process
        """
        assert exposure, "No exposure provided"

        image = exposure.getMaskedImage().getImage()

        ccd = afwCG.cast_Ccd(exposure.getDetector())

        for amp in ccd:
            if False:
                linear_threshold = amp.getElectronicParams().getLinearizationThreshold()
                linear_c = amp.getElectronicParams().getLinearizationCoefficient()
            else:
                linearizationCoefficient = self.config.linearizationCoefficient
                linearizationThreshold = self.config.linearizationThreshold

            if linearizationCoefficient == 0.0:     # nothing to do
                continue

            self.log.log(self.log.INFO,
                         "Applying linearity corrections to Ccd %s Amp %s" % (ccd.getId(), amp.getId()))

            if linearizationThreshold > 0:
                log10_thresh = math.log10(linearizationThreshold)

            ampImage = image.Factory(image, amp.getDataSec(), afwImage.LOCAL)

            width, height = ampImage.getDimensions()

            if linearizationThreshold <= 0:
                tmp = ampImage.Factory(ampImage, True)
                tmp.scaledMultiplies(linearizationCoefficient, ampImage)
                ampImage += tmp
            else:
                for y in range(height):
                    for x in range(width):
                        val = ampImage.get(x, y)
                        if val > linearizationThreshold:
                            val += val*linearizationCoefficient*(math.log10(val) - log10_thresh)
                            ampImage.set(x, y, val)
