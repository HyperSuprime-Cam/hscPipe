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


#    def illuminationCorrection(self, maskedImage, illumMaskedImage, illumscaling):
#
#        # common input test
#
#        maskedImage.scaledDivides(1./illumscaling, illumMaskedImage)

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
        # I am doing the measurement for QA separetely in the above lines, in this version.

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
        elif fittype == 'SUBARU': ## FH I will implement line-by-line subtraction
            offset = afwMath.makeStatistics(overscanData, afwMath.MEDIAN).getValue(afwMath.MEDIAN)
            maskedImage -= offset
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
        

