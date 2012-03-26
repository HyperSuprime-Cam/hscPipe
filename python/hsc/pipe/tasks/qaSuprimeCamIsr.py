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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import math, re
import numpy
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDetection
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg

## FH added for QA output
from lsst.ip.isr.isr import Isr
import lsst.afw.geom as afwGeom
#import hsc.onsite.measSeeingQa as qaSeeing
import hsc.onsite.qa.fitsthumb as qaFitsthumb


class QaSuprimeCamIsr(Isr):
    def __init__(self, display=False):
        self.display = display
    def createPsf(self, fwhm):
        """Make a PSF"""
        ksize = 4*int(fwhm) + 1
        return afwDetection.createPsf('DoubleGaussian', ksize, ksize, fwhm/(2*math.sqrt(2*math.log(2))))

    def calcEffectiveGain(self, maskedImage):
        im = afwImage.ImageF(maskedImage.getImage(), True)
        var = maskedImage.getVariance()
        im /= var
        medgain = afwMath.makeStatistics(im, afwMath.MEDIAN).getValue()
        meangain = afwMath.makeStatistics(im, afwMath.MEANCLIP).getValue()
        return medgain, meangain

    def calculateSdqaCcdRatings(self, maskedImage, metadata):
        metrics = {}
        metrics['nSaturatePix'] = 0
        metrics['nBadCalibPix'] = 0
        metrics['imageClipMean4Sig3Pass'] = None
        metrics['imageSigma'] = None
        metrics['imageMedian'] = None
        metrics['imageMin'] = None
        metrics['imageMax'] = None
        mask = maskedImage.getMask()
        badbitmask = mask.getPlaneBitMask('BAD')
        satbitmask = mask.getPlaneBitMask('SAT')
        intrpbitmask = mask.getPlaneBitMask('INTRP')
        sctrl = afwMath.StatisticsControl()
        sctrl.setNumIter(3)
        sctrl.setNumSigmaClip(4)
        sctrl.setAndMask(satbitmask | badbitmask | intrpbitmask)
        satmask = afwImage.MaskU(mask, True)
        badmask = afwImage.MaskU(mask, True)
        satmask &= satbitmask
        badmask &= badbitmask
        satmaskim = afwImage.ImageU(satmask.getBBox(afwImage.PARENT))
        satmaskim <<= satmask
        badmaskim = afwImage.ImageU(badmask.getBBox(afwImage.PARENT))
        badmaskim <<= badmask
        thresh = afwDetection.Threshold(0.5)
        fs = afwDetection.FootprintSet(satmaskim, thresh)
        for f in fs.getFootprints():
            metrics['nSaturatePix'] += f.getNpix()
        fs = afwDetection.FootprintSet(badmaskim, thresh)
        for f in fs.getFootprints():
            metrics['nBadCalibPix'] += f.getNpix()
        stats = afwMath.makeStatistics(maskedImage, afwMath.MEANCLIP | \
            afwMath.STDEVCLIP | afwMath.MEDIAN | afwMath.MIN |\
            afwMath.MAX, sctrl)
        metrics['imageClipMean4Sig3Pass'] = stats.getValue(afwMath.MEANCLIP)
        metrics['imageSigma'] = stats.getValue(afwMath.STDEVCLIP)
        metrics['imageMedian'] = stats.getValue(afwMath.MEDIAN)
        metrics['imageMin'] = stats.getValue(afwMath.MIN)
        metrics['imageMax'] = stats.getValue(afwMath.MAX)
        for k in metrics.keys():
            metadata.set(k, metrics[k])

    def calculateSdqaAmpRatings(self, maskedImage, metadata, biasBBox, dataBBox):
        metrics = {}
        metrics['nSaturatePix'] = 0
        metrics['overscanMean'] = None
        metrics['overscanStdDev'] = None
        metrics['overscanMedian'] = None
        metrics['overscanMin'] = None
        metrics['overscanMax'] = None
        trimmi = afwImage.MaskedImageF(maskedImage, dataBBox, False)
        biasmi = afwImage.MaskedImageF(maskedImage, biasBBox, False)
        mask = maskedImage.getMask()
        satbitmask = mask.getPlaneBitMask('SAT')
        sctrl = afwMath.StatisticsControl()
        sctrl.setAndMask(satbitmask)
        satmask = trimmi.getMask()
        satmask &= satbitmask
        satmaskim = afwImage.ImageU(satmask.getBBox(afwImage.PARENT))
        satmaskim <<= satmask
        thresh = afwDetection.Threshold(0.5)
        fs = afwDetection.FootprintSet(satmaskim, thresh)
        for f in fs.getFootprints():
            metrics['nSaturatePix'] += f.getNpix()
        stats = afwMath.makeStatistics(biasmi, afwMath.MEAN | \
            afwMath.STDEV | afwMath.MEDIAN | afwMath.MIN |\
            afwMath.MAX,sctrl)
        metrics['overscanMean'] = stats.getValue(afwMath.MEAN)
        metrics['overscanStdDev'] = stats.getValue(afwMath.STDEV)
        metrics['overscanMedian'] = stats.getValue(afwMath.MEDIAN)
        metrics['overscanMin'] = stats.getValue(afwMath.MIN)
        metrics['overscanMax'] = stats.getValue(afwMath.MAX)
        for k in metrics.keys():
            metadata.set(k, metrics[k])

    def interpolateDefectList(self, maskedImage, defectList, fwhm, fallbackValue=None):
        psf = self.createPsf(fwhm)
        if fallbackValue is None:
            fallbackValue = afwMath.makeStatistics(maskedImage.getImage(), afwMath.MEANCLIP).getValue()
        measAlg.interpolateOverDefects(maskedImage, psf, defectList, fallbackValue)
 
    def defectListFromFootprintList(self, fpList, growFootprints=1):
        defectList = measAlg.DefectListT()
        for fp in fpList:
            if growFootprints > 0:
                # if "True", growing requires a convolution
                # if "False", its faster
                fpGrow = afwDetection.growFootprint(fp, growFootprints, False)
            else:
                fpGrow = fp
            for bbox in afwDetection.footprintToBBoxList(fpGrow):
                defect = measAlg.Defect(bbox)
                defectList.push_back(defect)
        return defectList
   
    def maskPixelsFromDefectList(self, maskedImage, defectList, maskName='BAD'):
        # mask bad pixels
        mask = maskedImage.getMask()
        bitmask = mask.getPlaneBitMask(maskName)
        for defect in defectList:
            bbox = defect.getBBox()
            afwDetection.setMaskFromFootprint(mask, afwDetection.Footprint(bbox), bitmask)

    def getDefectListFromMask(self, maskedImage, maskName, growFootprints=1):
        mask = maskedImage.getMask()
        workmask = afwImage.MaskU(mask, True)
        workmask &= mask.getPlaneBitMask(maskName)
        thresh = afwDetection.Threshold(0.5)
        maskimg = afwImage.ImageU(workmask.getBBox(afwImage.PARENT))
        maskimg <<= workmask
        ds = afwDetection.FootprintSet(maskimg, thresh)
        fpList = ds.getFootprints()
        return self.defectListFromFootprintList(fpList, growFootprints)

    def makeThresholdMask(self, maskedImage, threshold, growFootprints=1, maskName = 'SAT'):
        if self.display:
            ds9.mtv(maskedImage, frame=0)

        # find saturated regions
        thresh = afwDetection.Threshold(threshold)
        ds = afwDetection.FootprintSet(maskedImage, thresh)
        fpList = ds.getFootprints()
        # set mask
        mask = maskedImage.getMask()
        bitmask = mask.getPlaneBitMask(maskName)
        if growFootprints > 0:
            for fp in fpList:
                fp = afwDetection.growFootprint(fp, growFootprints)
        afwDetection.setMaskFromFootprintList(mask, fpList, bitmask)

        return self.defectListFromFootprintList(fpList, growFootprints=0)

    def interpolateFromMask(self, maskedImage, fwhm, growFootprints = 1, maskName = 'SAT'):
        defectList = self.getDefectListFromMask(maskedImage, maskName, growFootprints)
        if 'INTRP' not in maskedImage.getMask().getMaskPlaneDict().keys():
            maskedImage.getMask.addMaskPlane('INTRP')
        psf = self.createPsf(fwhm)
        measAlg.interpolateOverDefects(maskedImage, psf, defectList)


    def saturationCorrection(self, maskedImage, saturation, fwhm, growFootprints=1, interpolate = True, maskName = 'SAT'):
        defectList = self.makeThresholdMask(maskedImage, saturation, grwoFootprints=growFootprints, maskName=maskName)
        if interpolate:
            measAlg.interpolateOverDefects(maskedImage, self.createPsf(fwhm), defectList)
        if self.display:
            ds9.mtv(maskedImage, frame=0)


    def biasCorrection(self, maskedImage, biasMaskedImage):

        maskedImage -= biasMaskedImage

    def darkCorrection(self, maskedImage, darkMaskedImage, expscaling, darkscaling):

        scale = expscaling / darkscaling
        maskedImage.scaledMinus(scale, darkMaskedImage)

    def updateVariance(self, maskedImage, gain):
        var = maskedImage.getVariance()
        var <<= maskedImage.getImage()
        var /= gain

#    def flatCorrection(self, maskedImage, flatMaskedImage, scalingtype, scaling = 1.0):
#        flatscaling = 1.0
#        # Figure out scaling from the data
#        # I'm not sure we should be doing this here, but maybe
#        if scalingtype == 'MEAN':
#            flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
#        elif scalingtype == 'MEDIAN':
#            flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
#        elif scalingtype == 'USER':
#            flatscaling = scaling
#        else:
#            raise pexExcept.LsstException, '%s : %s not implemented' % ("flatCorrection", scalingtype)
#        
#        maskedImage.scaledDivides(1./flatscaling, flatMaskedImage)
#
#        if self.display:
#            ds9.mtv(maskedImage, title="Flattened")

##== FH for QA output
    def flatCorrectionQa(self, maskedImage, flatMaskedImage, scalingtype, scaling = 1.0):
        flatscaling = 1.0
        # Figure out scaling from the data
        # I'm not sure we should be doing this here, but maybe
        if scalingtype == 'MEAN':
            flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEAN).getValue(afwMath.MEAN)
        elif scalingtype == 'MEDIAN':
            flatscaling = afwMath.makeStatistics(flatMaskedImage.getImage(), afwMath.MEDIAN).getValue(afwMath.MEDIAN)
        elif scalingtype == 'USER':
            flatscaling = scaling
        else:
            raise pexExcept.LsstException, '%s : %s not implemented' % ("flatCorrection", scalingtype)
        
        maskedImage.scaledDivides(1./flatscaling, flatMaskedImage)

        if self.display:
            ds9.mtv(maskedImage, title="Flattened")

        return
    
##== FH added for QA output
    def measureFlatnessImageQa(self, maskedImage, meshX=256, meshY=256, doClip=True, clipSigma=3, nIter=3):

        miSize = maskedImage.getDimensions()
        xmax = miSize[0] + int(meshX/2.) 
        ymax = miSize[1] + int(meshY/2.)
        nX = int(xmax / meshX)
        nY = int(ymax / meshY)
        skyLevels = numpy.zeros((nX,nY))

        # calcluating flatlevel over the subgrids 
        meshXHalf = int(meshX/2.)
        meshYHalf = int(meshY/2.)

        if doClip is True:
            sctrl = afwMath.StatisticsControl(clipSigma, nIter)

        for j in range(nY):
            yc = meshYHalf + j * meshY
            for i in range(nX):
                xc = meshXHalf + i * meshX
                xLLC = xc - meshXHalf
                yLLC = yc - meshYHalf
                xURC = xc + meshXHalf - 1
                yURC = yc + meshYHalf - 1
                
                bbox = afwGeom.Box2I(afwGeom.Point2I(xLLC, yLLC), afwGeom.Point2I(xURC, yURC))
                miMesh = maskedImage.Factory(maskedImage, bbox, afwImage.LOCAL)

                #clipSigma = 3.0
                #nIter = 3
                #stats = afwMath.makeStatistics(miMesh, afwMath.MEDIAN|afwMath.STDEVCLIP, sctrl)
                #stats = afwMath.makeStatistics(miMesh, afwMath.MEDIAN|afwMath.MEAN|afwMath.MEANCLIP, sctrl)
                stats = afwMath.makeStatistics(miMesh, afwMath.MEDIAN|afwMath.MEAN|afwMath.MEANCLIP, sctrl)                

                if doClip is True:
                    skyLevels[i, j] = stats.getValue(afwMath.MEANCLIP)
                else:
                    skyLevels[i, j] = stats.getValue(afwMath.MEAN)
                #skyLevels[i, j] = stats.getValue(afwMath.MEDIAN)
                #skySigma[i, j] = stats.getValue(afwMath.STDEVCLIP)                

        skyMedian = numpy.median(skyLevels)
        flatness =  (skyLevels - skyMedian) / skyMedian        
        flatness_rms = numpy.std(flatness)
        flatness_min = flatness.min()
        flatness_max = flatness.max() 
        flatness_pp = flatness_max - flatness_min

        return (flatness, flatness_pp, flatness_min, flatness_max, flatness_rms, skyMedian, nX, nY)


    def illuminationCorrection(self, maskedImage, illumMaskedImage, illumscaling):

        # common input test

        maskedImage.scaledDivides(1./illumscaling, illumMaskedImage)



    def trimAmp(self, exposure, trimBbox=None):
        """
        This returns a new Exposure that is a subsection of the input exposure.

        NOTE : do we need to deal with the WCS in any way, shape, or form?
        """
        if trimBbox is not None:
            return exposureFactory(exposure, trimBbox, LOCAL)
        else:
            amp = cameraGeom.cast_Amp(exposure.getDetector())
            return exposureFactory(exposure, amp.getDiskDataSec(false), LOCAL)
        # n.b. what other changes are needed here?
        # e.g. wcs info, overscan, etc

#    def overscanCorrection(self, maskedImage, overscanData, fittype='MEDIAN', polyorder=1, imageFactory=afwImage.ImageF):
#        """
#        """
#        typemap = {afwImage.ImageU:numpy.uint16, afwImage.ImageI:numpy.int32, afwImage.ImageF:numpy.float32, afwImage.ImageD:numpy.float64}
#
#        # what type of overscan modeling?
#        offset = 0
#        if fittype == 'MEAN':
#            offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
#            maskedImage -= offset
#        elif fittype == 'MEDIAN':
#            offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
#            maskedImage -= offset
#        elif fittype == 'POLY':
#            biasArray = overscanData.getArray()
#            #Assume we want to fit along the long axis
#            aind = numpy.argmin(biasArray.shape)
#            find = numpy.argmin(biasArray.shape)
#            fitarr = numpy.median(biasArray, axis=aind)
#            coeffs = numpy.polyfit(range(len(fitarr)), fitarr, deg=polyorder)
#            offsets = numpy.polyval(coeffs, range(len(fitarr)))
#            width, height = maskedImage.getDimensions()
#            offarr = numpy.zeros((height, width), dtype = typemap[imageFactory])
#            if aind == 1:
#                for i in range(len(offsets)):
#                    offarr[i] = offsets[i]
#            elif aind == 0:
#                offarr = offarr.T
#                for i in range(len(offsets)):
#                    offarr[i] = offsets[i]
#                offarr = offarr.T
#            else:
#                raise pexExcept.LsstException, "Non-2D array returned from MaskedImage.getArray()"
#            im = afwImage.makeImageFromArray(offarr)
#            maskedImage -= im 
#        else:
#            raise pexExcept.LsstException, '%s : %s an invalid overscan type' % ("overscanCorrection", fittype)

    ##== FH added for QA output
    def overscanCorrectionQa(self, maskedImage, overscanData, fittype='MEDIAN', polyorder=1, imageFactory=afwImage.ImageF):
        """
        """
        typemap = {afwImage.ImageU:numpy.uint16, afwImage.ImageI:numpy.int32, afwImage.ImageF:numpy.float32, afwImage.ImageD:numpy.float64}

        # measurement of overscan levels for QA, using median with outlier clipping
        clipSigma = 3.0
        nIter = 3
        sctrl = afwMath.StatisticsControl(clipSigma, nIter)
        stats = afwMath.makeStatistics(overscanData, afwMath.MEDIAN|afwMath.STDEVCLIP, sctrl)
        osLevel = stats.getValue(afwMath.MEDIAN)
        osSigma = stats.getValue(afwMath.STDEVCLIP)

        #print "QA overscan: osLevel: %f" % osLevel
        #print "QA overscan: osSigma: %f" % osSigma

        # n.b. by fh: below part performs measurement of overscan level again 
        # separately from the above measurement for QA, which is a duplicated operation.
        # But, to use a single algorithm consistently for all frames for QA,
        # I am doing the measurement for QA separetely in the above, in this version.

        # what type of overscan modeling?
        offset = 0
        if fittype == 'MEAN':
            offset = afwMath.makeStatistics(overscanData, afwMath.MEAN).getValue(afwMath.MEAN)
            maskedImage -= offset
        elif fittype == 'MEDIAN':
            offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
            maskedImage -= offset
        elif fittype == 'POLY':
            biasArray = overscanData.getArray()
            #Assume we want to fit along the long axis
            aind = numpy.argmin(biasArray.shape)
            find = numpy.argmin(biasArray.shape)
            fitarr = numpy.median(biasArray, axis=aind)
            coeffs = numpy.polyfit(range(len(fitarr)), fitarr, deg=polyorder)
            offsets = numpy.polyval(coeffs, range(len(fitarr)))
            width, height = maskedImage.getDimensions()
            offarr = numpy.zeros((height, width), dtype = typemap[imageFactory])
            if aind == 1:
                for i in range(len(offsets)):
                    offarr[i] = offsets[i]
            elif aind == 0:
                offarr = offarr.T
                for i in range(len(offsets)):
                    offarr[i] = offsets[i]
                offarr = offarr.T
            else:
                raise pexExcept.LsstException, "Non-2D array returned from MaskedImage.getArray()"

            im = afwImage.makeImageFromArray(offarr)
            maskedImage -= im

        else:
            raise pexExcept.LsstException, '%s : %s an invalid overscan type' % ("overscanCorrection", fittype)

        return osLevel, osSigma

    ##== FH added for QA output
    def writeFitsImageQa(self, exposure, outfile):
        """
        writing out exposure to an FITS image named outfile. 
        """
        exposure.getMaskedImage().writeFits(outfile)

    ##== FH added for QA output
    def writeSnapshotImageQa(self, exposure, outfile, format='png', width=500, height=0):
        """
        writing out exposure to a snapshot file named outfile in the given image format and size.  
        """
        qaFitsthumb.createFitsThumb(exposure.getMaskedImage().getImage(), outfile, format, width, height, True)
        

    def fringeCorrection(self, maskedImage, fringe):

        raise pexExcept.LsstException, '%s not implemented' % ("ipIsr.fringCorrection")


    def pupilCorrection(self, maskedImage, pupil):

        raise pexExcept.LsstException, '%s not implemented' % (stageName)
'''
class Linearization(object):
    def __init__(self, linearityFile=None):
        if linearityFile is not None:
            self.readFile(linearityFile)
        else:
            self.makeLinearReplace()
    def getImageFactoryFromExposure(self, exposure):
        return exposure.getMaskedImage().getImage().Factory
    def apply(self, exposure):
        self.getImageFactoryFromExposure(exposure)
        if type is "LUT":
            imageData = exposure.getMaskedImage(). 
        mi = exposure.getMaskedImage()
        	
    def readFile(self, filename):
    def writeFile(self, filename):
    def makeLinearReplace(self):
        self.type = 
    def makeLinearMult(self):
'''
