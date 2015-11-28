#!/usr/bin/env python
"""
TODO:
 This module is to select stars (PSF-like sources) based on a sequence
 which comprises PSF sources seen in the 'mag vs fwhm' 2-d space.
  Application of the same strategy in seeing measurement in the SC-RCM running
  for Suprime-Cam observation (publised in Furusawa+2011,PASJ,63,S581).
  The advantage is to determine a decent limiting magnitude range for
   extraction of PSF sequence in the 'mag vs fwhm' space, to avoid
  contaminations from too many faint galaxies or cosmic rays.
  Initially implemented for onsite QA output, and now being ported for
  general use.
"""

import os, os.path
import math
import lsst.pex.config as pexConfig
import lsst.pex.logging as pexLog
import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as geomEllip
import lsst.meas.algorithms as measAlg

from lsst.pipe.base import Task, Struct
import lsst.afw.display.ds9 as ds9

import numpy
import errno

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#plt.switch_backend('Agg')

import matplotlib.figure as figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
from matplotlib import patches as patches


class SizeMagnitudeMitakaStarSelectorConfig(pexConfig.Config):
    fracSrcIni = pexConfig.Field(
        dtype = float,
        doc = 'What fraction of sources from the brightest is to be included for initial guess of seeing to avoid cosmic rays which dominate faint magnitudes', 
        default = 0.15, # is good for SC with 5-sigma detection. might be a bit too large in blue band
        )
    fwhmMin  = pexConfig.Field(
        dtype = float,
        doc = 'Minimum fwhm allowed in estimation of seeing (pix)',
        default = 1.5,
        )
    fwhmMax  = pexConfig.Field(
        dtype = float,
        doc = 'Maxmum fwhm allowed in estimation of seeing (pix)',
        default = 12.0,
        )
    nbinMagHist = pexConfig.Field(
        dtype = int,
        doc = 'Number of bins for number counting as a fn of instrumnetal mag',
        default = 80,
        )
    magMinHist = pexConfig.Field(
        dtype = float,
        doc = 'Brightest mag for number counting as a fn of instrumnetal mag',
        default = -20.0,
        )
    magMaxHist = pexConfig.Field(
        dtype = float,
        doc = 'Faintest mag for number counting as a fn of instrumnetal mag',
        default = 0.0,
        )
    statAlgRoughFwhm = pexConfig.ChoiceField(
        doc = 'Statistical algorithm to derive rough Fwhm in the 1st step seeing estimation',
        dtype = str,  default = "MEDIAN",
        allowed = {
        "MEDIAN" : "median of sample",
        "MEANCLIP" : "clipped mean of sample with 3-sigma clip + 3-times iteration",
        }
        )
    nBrightSampleRoughFwhm = pexConfig.Field(
        dtype = int,
        doc = 'Number of brightest (non-saturated) objects which are used to determine rough-interim seeing',
        #default = 50, # seems to be too small to capture lager fwhm sources?
        default = 30,
        )
    nSmallSampleRoughFwhm = pexConfig.Field(
        dtype = int,
        doc = 'Number of smallest objects which are used to determine rough-interim seeing',
        #default = 30, # seems to be too small to capture lager fwhm sources?
        default = 50,
        #default = 100,
        )
    fwhmMarginFinal = pexConfig.Field(
        dtype = float,
        doc = 'How many pixels around the peak are used for calculating scatter of psf candidates',
        #default = 0.5, # --> too tight to assess width of psf sequence?
        default = 0.75, # seems to be ok for SC CCD=0,a bit tight =8
        #default = 1.0, # seems to be ok for SC CCD=0,a bit tight =8
        #default = 1.5,
        )
    fwhmMarginNsigma = pexConfig.Field(
        dtype = float,
        doc = 'How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence',
        #default = 1.0,
        #default = 1.2, # good but not optimal for broad psf sequence
        #default = 1.5, 
        #default = 2.0, # a bit broad for some data
        ###default = 2.5, # a bit broad
        default = 3., # broad
        )
    psfSeqStatNsigma = pexConfig.Field(
        dtype = float,
        doc = 'How many sigmas around the peak fwhm are used for calculating statistics of PSF sequence',
        default = 3.,
        )
    psfSeqStatNiter = pexConfig.Field(
        dtype = int, 
        doc = 'How many times do we iterate calculating statistics of PSF sequence',
        default = 3,
        )
    magLimitFaintExtension = pexConfig.Field(
        dtype = float,
        doc = 'How many magnitudes to extend the faint-end limit for extracting PSF sources, from the base magnitude determined by fracSrcIni.',
        default = 0.0, # 0.0 would make clean sample
        )
    fwhmBinSize = pexConfig.Field(
        dtype = float,
        doc = "Bin size of FWHM histogram",
        default = 0.2,
        )
    doPlots = pexConfig.Field(
        dtype = bool,
        doc = "Make plots?",
        default = True
        )
    gridSize = pexConfig.Field(
        dtype = float,
        doc = "Size of grid (pixels)",
        default = 1024
        )
    kernelSize = pexConfig.Field(
        doc = "size of the kernel to create",
        dtype = int,
        default = 21,
    )
    borderWidth = pexConfig.Field(
        doc = "number of pixels to ignore around the edge of PSF candidate postage stamps",
        dtype = int,
        default = 0,
    )
    doUndistort = pexConfig.Field(
        dtype = bool,
        doc = "Undistort when evaluating the 2nd moments of sources?",
        default = False,
    )

class SizeMagnitudeMitakaStarSelector(object):
    """
    """
    ConfigClass = SizeMagnitudeMitakaStarSelectorConfig

    def __init__(self, config, schema=None, metadata=None, **kwrgs):
        """
        Construct a star selector that uses second moments and instrumental magnitude
        This is a naive algorithm and should be used with caution.

        @param[in] config: An instance of SizeMagnitudeMitakaStarSelectorConfig
        @param[in,out] schema: An afw.table.Schema to register the selector's flag field.
        If None, the sources will not be modified.
        """
        if not config:
            config = SizeMagnitudeMitakaStarSelector.ConfigClass()
        self.config = config

        self._kernelSize  = config.kernelSize
        self._borderWidth = config.borderWidth

        if schema is not None:
            starflagName = "classification.mitakastar"
            if starflagName in schema.getNames():
                self._key = schema.find("classification.mitakastar").key
            else:
                self._key = schema.addField("classification.mitakastar", type="Flag",
                                            doc="selected as a star by sizeMagnitudeMitakaStarSelector")

        #self.config.doPlots = False
        self.debugFlag = False

        # In the case of Task-based class, these two are created in the parent class' __init__
        self.log = pexLog.Log.getDefaultLog()
        if metadata is not None:
            self.metadata = metadata
        else:
            self.metadata = dafBase.PropertySet()


    def selectStars(self, exposure, catalog, dataRef=None, outputStruct=False):
        """
        Return a list of PSF candidates that represent likely stars

        A list of PSF candidates may be used by a PSF fitter to construct a PSF

        @param[in] dataRef: dataRef passed by the pipeline including reference to a butler instance
        @param[in] exposure: the exposure containing the sources
        @param[in] catalog: a source list containing sources that may be stars
        @return psfCandidateList: a list of PSF candidates.
        """

        self.log.info("Mitaka StarSelector has been called.")

        import lsstDebug
        display = lsstDebug.Info(__name__).display
        displayExposure = lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells
        pauseAtEnd = lsstDebug.Info(__name__).pauseAtEnd               # pause when done

        #display=True
        #displayExposure= True

        mi = exposure.getMaskedImage()

        if display:
            frames = {}
            if displayExposure:
                frames["displayExposure"] = 1
                ds9.mtv(mi, frame=frames["displayExposure"], title="PSF candidates")

        # screening sources
        goodData = self.getGoodSources(catalog, exposure)

        if display and displayExposure:
            with ds9.Buffering():
                for i in range(len(goodData.xListAll)):
                    if i in goodData.indicesSourcesFwhmRange:
                        ctype = ds9.GREEN # good
                    elif i in goodData.indicesGoodSources:
                        ctype = ds9.MAGENTA # rejected by fwhm search range
                    else:
                        ctype = ds9.RED # bad; saturation etc.
                    xc = goodData.xListAll[i]
                    yc = goodData.yListAll[i]
                    ds9.dot("o", xc, yc, frame=frames["displayExposure"], ctype=ctype)

        # determine mag limit for clean psf candidates, dynamically based on the image
        magLimPsfSeq = self.getMagLimit(dataRef, goodData, exposure)
        if magLimPsfSeq is None:
            if outputStruct:
                return None, None
            else:
                return None

        # getting a decent median fwhm
        #   fwhm = -9999.0 is returned in cases of failure 
        fwhmRough = self.getFwhmRough(dataRef, goodData, magLimPsfSeq, exposure)

        # extract a list of psf-like sources
        # goodData is updated as same as dataPsfLike in getStarCandidateList()
        if self.config.doUndistort:
            fwhmRough = goodData.fwhmUndistRough
        dataPsfLike = self.getStarCandidateList(dataRef, goodData, fwhmRough, magLimPsfSeq)

        if dataPsfLike is None:
            self.log.warn("psf-like candadate list from seeing rough estimation is empty")
            if outputStruct:
                return None, None
            else:
                return None

        # merging local metadata into exposure metadata
        # n.g., when called from QaTask, this is done again in the QaTask.run() in current version.
        #exposure.getMetadata().combine(self.metadata)
        metadata = exposure.getMetadata()
        for key in self.metadata.names():
            metadata.set(key, self.metadata.get(key))
            #print '*** local.metadata: %s = %s is set to exposure' % (key, str(self.metadata.get(key)))

        if len(catalog) != len(dataPsfLike.xListAll):
            self.log.warn("Number of catalog sources %d mismatch with gooddata list %d" % (len(catalog), len(dataPsfLike.xListAll)))
            if outputStruct:
                return None, None
            else:
                return None

        xCcdSize, yCcdSize = exposure.getWidth(), exposure.getHeight()
        psfCandidateList = []

        with ds9.Buffering():
            for i in dataPsfLike.indicesSourcesPsfLikeRobust:
                source = catalog[i]
                try:
                    psfCandidate = measAlg.makePsfCandidate(source, exposure)
                    # The setXXX methods are class static, but it's convenient to call them on
                    # an instance as we don't know Exposure's pixel type
                    # (and hence psfCandidate's exact type)
                    if psfCandidate.getWidth() == 0:
                        psfCandidate.setBorderWidth(self._borderWidth)
                        psfCandidate.setWidth(self._kernelSize + 2*self._borderWidth)
                        psfCandidate.setHeight(self._kernelSize + 2*self._borderWidth)

                    # to check if the psf image sits within the CCD pixel area
                    xPsf, yPsf = psfCandidate.getXCenter(), psfCandidate.getYCenter()
                    xPsfSize, yPsfSize = psfCandidate.getWidth(), psfCandidate.getHeight()
                    if not (xPsfSize*0.5 < xPsf and xPsf <= (xCcdSize-xPsfSize*0.5) and \
                            yPsfSize*0.5 < yPsf and yPsf <= (yCcdSize-yPsfSize*0.5)):
                        continue

                    im = psfCandidate.getMaskedImage().getImage()
                    maxVal = afwMath.makeStatistics(im, afwMath.MAX).getValue()
                    if not numpy.isfinite(maxVal):
                        continue
                    if self._key is not None:
                        source.set(self._key, True)
                    psfCandidateList.append(psfCandidate)

                    symb, ctype = "+", ds9.GREEN
                except Exception as err:
                    symb, ctype = "o", ds9.RED
                    print "FH", err
                    pass # FIXME: should log this!

                if display and displayExposure:
                    ds9.dot(symb, source.getX() - mi.getX0(), source.getY() - mi.getY0(),
                            size=4, frame=frames["displayExposure"], ctype=ctype)

        if display and pauseAtEnd:
            raw_input("Continue? y[es] p[db] ")

        if outputStruct:
            return psfCandidateList, dataPsfLike
        else:
            return psfCandidateList


    def isGoodSource(self, source, keySaturationCenterFlag, keySaturationFlag):
        Ixx = source.getIxx()
        Iyy = source.getIyy()
        Ixy = source.getIxy()
        sigma = math.sqrt(0.5*(Ixx+Iyy))
        fwhm = 2.*math.sqrt(2.*math.log(2.)) * sigma # (pix) assuming Gaussian
        saturationFlag = source.get(keySaturationFlag)
        saturationCenterFlag = source.get(keySaturationCenterFlag)
        isSaturated = (saturationFlag | saturationCenterFlag)
        if self.debugFlag:
            print ('objFlag SatAny %x  SatCenter %x  isSaturated %x' %
                   (saturationFlag, saturationCenterFlag, isSaturated))

        fluxAperture = source.getApFlux()
        fluxGauss = source.getModelFlux()
        fluxPsf = source.getPsfFlux()
        fluxForSeeing = fluxAperture

        if  isSaturated or math.isnan(fwhm) or math.isnan(Ixx) or math.isnan(fluxForSeeing) or \
           fwhm <= 0 or fluxForSeeing <= 0 or (Ixx*Iyy*Ixy == 0):
            return False
        else:
            return True

    def getEllipticityFromSecondmoments(self, Ixx, Iyy, Ixy):
        VerySmallValue = 1.0e-10
        if True: # definition by SExtractor
            Val1 = 0.5*(Ixx+Iyy)
            Ixx_Iyy = Ixx-Iyy
            Val2 = 0.25*Ixx_Iyy*Ixx_Iyy + Ixy*Ixy
            if Val2 >= 0 and (Val1-math.sqrt(Val2)) > 0:
                aa = math.sqrt( Val1 + math.sqrt(Val2) )
                bb = math.sqrt( Val1 - math.sqrt(Val2) )
                ell =  1. - bb/aa
                if math.fabs(Ixx_Iyy) > VerySmallValue:
                    tanVal = 2.0 * Ixy / Ixx_Iyy
                    if tanVal >= 0:
                        if Ixx_Iyy > 0: # elongation toward x
                            ellPa = 0.5*math.atan(tanVal)
                        else: # elongation toward y
                            ellPa = 0.5*(math.atan(tanVal)+math.pi)
                    else:
                        if Ixx_Iyy > 0: # elongation toward x
                            ellPa = 0.5*(math.atan(tanVal)+2.0*math.pi)
                        else: # elongation toward y
                            ellPa = 0.5*(math.atan(tanVal)+math.pi)
                else: # source is too round to estimate PA of elongation
                    ellPa = -9999.0
            else:
                ell = -9999.0
                ellPa = -9999.0
                aa = -9999.0
                bb = -9999.0
        if True: # definition by Kaiser
            # e=sqrt(e1^2+e2^2) where e1=(Ixx-Iyy)/(Ixx+Iyy), e2=2Ixy/(Ixx+Iy)
            # SExtractor's B/A=sqrt((1-e)/(1+e)), ell=1-B/A
            fabs_Ixx_p_Iyy = math.fabs(Ixx+Iyy)
            if fabs_Ixx_p_Iyy > VerySmallValue:
                e1 = (Ixx-Iyy)/(Ixx+Iyy)
                ### if e1 > 0: ### XXX need to confirm this condition is unnecessary
                e2 = 2.0*Ixy/(Ixx+Iyy)
                ell_e1e2 = math.sqrt(e1*e1 + e2*e2)
                fabs_Ixx_m_Iyy = math.fabs(Ixx-Iyy)
                if fabs_Ixx_m_Iyy > VerySmallValue:
                    tanVal = 2.0 * Ixy / Ixx_Iyy
                    if tanVal >= 0:
                        if Ixx_Iyy > 0: # elongation toward x
                            ellPa_e1e2 = 0.5*math.atan(tanVal)
                        else: # elongation toward y
                            ellPa_e1e2 = 0.5*(math.atan(tanVal)+math.pi)
                    else:
                        if Ixx_Iyy > 0: # elongation toward x
                            ellPa_e1e2 = 0.5*(math.atan(tanVal)+2.0*math.pi)
                        else: # elongation toward y
                            ellPa_e1e2 = 0.5*(math.atan(tanVal)+math.pi)
                else: # source is too round to estimate PA of elongation
                    ellPa_e1e2 = -9999.0
            else:
                e1 = -9999.0
                e2 = -9999.0
                ell_e1e2 = -9999
                ellPa_e1e2 = -9999.0


        if -90.0 <= ellPa and ellPa <= 90.0:
            #ellPa = 90. - ellPa ## definition of PA to be confirmed
            ellPa = math.degrees(ellPa)
            pass

        return Struct(
            ell = ell,
            aa = aa,
            bb = bb,
            ellPa = ellPa,
            e1 = e1,
            e2 = e2,
            ell_e1e2 = ell_e1e2,
            )

    def getGoodSources(self, catalog, exposure):
        """screening sources to get not-dirty sources from catalog"""

        # -- Filtering sources with rough min-max fwhm
        xListAll = []
        yListAll = []
        magListAll = []
        fwhmListAll = []
        IxxListAll = []
        IyyListAll = []
        IxyListAll = []
        fwhmUndistListAll = []
        IxxUndistListAll = []
        IyyUndistListAll = []
        IxyUndistListAll = []

        indicesGoodSources = [] # indices of sources which are not dirty
        indicesSourcesFwhmRange = [] # indices of sources in acceptable fwhm range

        # below not relevant to star selector but required for QA
        ellListAll = []
        ellPaListAll = []
        AEllListAll = []
        BEllListAll = []

        e1ListAll = []
        e2ListAll = []
        elle1e2ListAll = []

        # Undistorting moments when requested
        if self.config.doUndistort:
            detector = exposure.getDetector()
            distorter = None
            xy0 = afwGeom.Point2D(0,0)
            if not detector is None:
                # Note: we use getCenter() instead of getCenterPixel() because getCenterPixel() assumes
                # that CCDs are laid out in a regular grid, which may not be true (e.g., HSC).
                pixSize = detector.getPixelSize()
                cPix = detector.getCenter().getPixels(pixSize)
                detSize = detector.getSize().getPixels(pixSize)
                xy0.setX(cPix[0] - int(0.5*detSize[0]))
                xy0.setY(cPix[1] - int(0.5*detSize[1]))
                distorter = detector.getDistortion()

        catalogSchema = catalog.table.getSchema()
        nameSaturationFlag = 'flags.pixel.saturated.any'
        nameSaturationCenterFlag = 'flags.pixel.saturated.center'
        keySaturationFlag = catalogSchema.find(nameSaturationFlag).key
        keySaturationCenterFlag = catalogSchema.find(nameSaturationCenterFlag).key

        for iseq, source in enumerate(catalog):

            xc = source.getX()
            yc = source.getY()
            Ixx = source.getIxx() # luminosity-weighted 2nd moment of pixels
            Iyy = source.getIyy()
            Ixy = source.getIxy()
            sigma = math.sqrt(0.5*(Ixx+Iyy))
            fwhm = 2.*math.sqrt(2.*math.log(2.)) * sigma # (pix) assuming Gaussian
            fwhmToEval = fwhm

            fluxAper = source.getApFlux()
            fluxErrAper =  source.getApFluxErr()
            fluxGauss = source.getModelFlux()
            fluxErrGauss =  source.getModelFluxErr()
            fluxPsf = source.getPsfFlux()
            fluxErrPsf = source.getPsfFluxErr()

            # which flux do you use:
            fluxForSeeing = fluxAper
            mag = -2.5*numpy.log10(fluxForSeeing)

            xListAll.append(xc)
            yListAll.append(yc)
            IxxListAll.append(Ixx)
            IyyListAll.append(Iyy)
            IxyListAll.append(Ixy)
            magListAll.append(mag)
            fwhmListAll.append(fwhm)

            if self.config.doUndistort:
                xpix, ypix = xc + xy0.getX(), yc + xy0.getY()
                p = afwGeom.Point2D(xpix, ypix)
                m = distorter.undistort(p, geomEllip.Quadrupole(Ixx, Iyy, Ixy), detector)
                IxxUndist, IyyUndist, IxyUndist = m.getIxx(), m.getIyy(), m.getIxy()
                sigma = math.sqrt(0.5*(IxxUndist+IyyUndist))
                fwhmUndist = 2.*math.sqrt(2.*math.log(2.)) * sigma
                fwhmToEval = fwhmUndist

                IxxUndistListAll.append(IxxUndist)
                IyyUndistListAll.append(IyyUndist)
                IxyUndistListAll.append(IxyUndist)
                fwhmUndistListAll.append(fwhmUndist)

            # these for QA
            ellRet = self.getEllipticityFromSecondmoments(Ixx, Iyy, Ixy)
            ellListAll.append(ellRet.ell)
            ellPaListAll.append(ellRet.ellPa)
            AEllListAll.append(ellRet.aa)
            BEllListAll.append(ellRet.bb)
            e1ListAll.append(ellRet.e1)
            e2ListAll.append(ellRet.e2)
            elle1e2ListAll.append(ellRet.ell_e1e2)

            # checking if this source has valid measurements
            if not self.isGoodSource(source, keySaturationCenterFlag, keySaturationFlag):
                continue
            else:
                indicesGoodSources.append(iseq)

                # checking if this source is sitting within the search range of Fwhm?
                if self.config.fwhmMin < fwhmToEval and fwhmToEval < self.config.fwhmMax:
                    indicesSourcesFwhmRange.append(iseq)
                else:
                    continue

        return Struct(
            magListAll = numpy.array(magListAll),
            fwhmListAll = numpy.array(fwhmListAll),
            xListAll = numpy.array(xListAll),
            yListAll = numpy.array(yListAll),
            IxxListAll = numpy.array(IxxListAll),
            IyyListAll = numpy.array(IyyListAll),
            IxyListAll = numpy.array(IxyListAll),
            ellListAll = numpy.array(ellListAll),
            ellPaListAll = numpy.array(ellPaListAll),
            AEllListAll = numpy.array(AEllListAll),
            BEllListAll = numpy.array(BEllListAll),
            e1ListAll = numpy.array(e1ListAll),
            e2ListAll = numpy.array(e2ListAll),
            elle1e2ListAll = numpy.array(elle1e2ListAll),
            indicesGoodSources = indicesGoodSources,
            indicesSourcesFwhmRange = indicesSourcesFwhmRange,
            IxxUndistListAll = numpy.array(IxxUndistListAll),
            IyyUndistListAll = numpy.array(IyyUndistListAll),
            IxyUndistListAll = numpy.array(IxyUndistListAll),
            fwhmUndistListAll = numpy.array(fwhmUndistListAll),
            )

    def getMagLimit(self, dataRef, data, exposure):
        """
        Derive the normalized cumulative magnitude histogram and estimate good limit mag
        for extracting psf-like candidates based on the histogram
        """

        #import pudb; pudb.set_trace()

        magListFwhmRange = data.magListAll[data.indicesSourcesFwhmRange]
        # n.b. below magHist[1] is sorted so index=0 is brightest
        magHist = numpy.histogram(magListFwhmRange, range=(self.config.magMinHist, self.config.magMaxHist),
                                  bins=self.config.nbinMagHist)

        sumAll = magHist[0].sum()
        #print '*** sumAll: ', sumAll
        self.log.info("QaSeeing: total number of objects in the first dataset (sumAll): %d" % sumAll)

        magCumHist = list(magHist)
        magCumHist[0] = (magCumHist[0].cumsum())
        #print '------ mag cum hist ------'
        #print magCumHist
        magCumHist[0] = magCumHist[0] / float(sumAll)
        #print '------ normalized mag cum hist ------'
        #print magCumHist

        # -- Estimating mag limit based no the cumlative mag histogram
        magLim = None
        zpFrame = None

        # handling slight dependency on filter kinds
        filterName = exposure.getFilter().getName()
        if filterName in ['g', 'B']:
            fracSrcIni = self.config.fracSrcIni * (2./3.)
        else:
            fracSrcIni = self.config.fracSrcIni
        self.log.info("Filter: %s -> fraction for guessing magnitude limit is set to: %f" % (filterName, fracSrcIni))

        for i, cumFraction in enumerate(magCumHist[0]):
            if cumFraction >= fracSrcIni:
                magLim = magCumHist[1][i] # magLim is the mag which exceeds the cumulative n(m) of 0.15
                break
        if magLim is None or numpy.isnan(magLim):
            magLim = 99.0
            magLimCalib = 99.0
            fluxLim = -9999.0
            self.log.log(self.log.WARN, "Error: cumulative magnitude histogram does not exceed fracSrcIni: %f" % fracSrcIni)
        else:
            fluxLim = numpy.power(10, -0.4*magLim)
            try:
                zpFrame = exposure.getCalib().getMagnitude(1.0) # zp(mag/ADU/exptime)
            except: # photocal yet to be done or failure in photocal
                zpFrame = None
            if zpFrame is None or numpy.isnan(zpFrame):
                magLimCalib = 99.0
            else:
                magLimCalib = magLim + zpFrame

        #print '*** magLim: ', magLim
        self.log.info("Mag limit auto-determined: %5.2f (%5.2f by calib) or %f (ADU)" % (magLim, magLimCalib, fluxLim))
        self.metadata.set("magLim", magLim)
        self.metadata.set("magLimCalib", magLimCalib)
        self.metadata.set("fluxLim", fluxLim)

        if self.debugFlag:
            self.log.logdebug("QaSeeing: magHist: %s" % magHist)
            self.log.logdebug("QaSeeing: cummurative magHist: %s" % magCumHist)

        if self.config.doPlots and dataRef is not None:
            fig = figure.Figure()
            canvas = FigCanvas(fig)
            #fig = plt.figure()
            pltMagHist = fig.add_subplot(2,1,1)
            pltMagHist.hist(magListFwhmRange, bins=magHist[1], orientation='vertical')
            pltMagHist.set_title('differential / cumulative histograms of magnitudes')
            #pltMagHist.set_xlabel('magnitude instrumental')
            pltMagHist.set_xticklabels([])
            pltMagHist.set_ylabel('number of samples')
            pltMagHist.legend()

            pltCumHist = fig.add_subplot(2,1,2)
            # histtype=bar,barstacked,step,stepfilled
            pltCumHist.hist(magListFwhmRange, bins=magCumHist[1], normed=True, cumulative=True,
                            orientation='vertical', histtype='step')
            xx = [magLim, magLim]; yy = [0, 1]
            pltCumHist.plot(xx, yy, linestyle='dashed', label='mag limit') # solid,dashed,dashdot,dotted

            if zpFrame and not numpy.isnan(zpFrame):
                pltCumHist2 = pltCumHist.twiny() # another x axis with the calibrated mag
                xlimCalib = [ x + zpFrame for x in pltCumHist.get_xlim() ]
                pltCumHist2.set_xlim(xlimCalib)
                pltCumHist2.set_ylim(pltCumHist.get_ylim())
                pltCumHist2.plot([-100,100],[-1,-1])

            pltCumHist.set_xlabel('magnitude (bottom:isnt / top:calib)')
            pltCumHist.set_ylabel('Nsample scaled to unity')
            pltCumHist.legend()

            fname = getFilename(dataRef, "plotMagHist")
            fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                        format='png', transparent=False, bbox_inches=None, pad_inches=0.1)
            del fig
            del pltMagHist
            del pltCumHist
            del canvas
            if zpFrame and not numpy.isnan(zpFrame):
                del pltCumHist2

        return magLim

    def getFwhmRough(self, dataRef, data, magLimSeq, exposure):
#    def getFwhmRough(self, magListAll, fwhmListAll, indicesSourcesFwhmRange, magLim):

        """Estimating roughly-estimated FWHM for sources with mag < magLim"""
        # saturated sources are already filtered in data.indicesSourcesFwhmRange
        indicesSourcesForRoughFwhm = [i for i in data.indicesSourcesFwhmRange if data.magListAll[i] < magLimSeq ]

        magListForRoughFwhm = data.magListAll[indicesSourcesForRoughFwhm]
        fwhmListForRoughFwhm = data.fwhmListAll[indicesSourcesForRoughFwhm]
        if self.config.doUndistort:
            fwhmUndistListForRoughFwhm = data.fwhmUndistListAll[indicesSourcesForRoughFwhm]
        else:
            fwhmUndistListForRoughFwhm = []

        data.indicesSourcesForRoughFwhm = indicesSourcesForRoughFwhm
        data.magListForRoughFwhm = magListForRoughFwhm
        data.fwhmListForRoughFwhm = fwhmListForRoughFwhm
        data.fwhmUndistListForRoughFwhm = fwhmUndistListForRoughFwhm

        # below for the QA plots
        #xListForRoughFwhm = data.xListAll[indicesSourcesForRoughFwhm]
        #yListForRoughFwhm = data.yListAll[indicesSourcesForRoughFwhm]
        #IxxListForRoughFwhm = data.IxxListAll[indicesSourcesForRoughFwhm]
        #IyyListForRoughFwhm = data.IyyListAll[indicesSourcesForRoughFwhm]
        #IxyListForRoughFwhm = data.IxyListAll[indicesSourcesForRoughFwhm]

        self.log.logdebug("nBritestSampleRoughFwhm: %d" % self.config.nBrightSampleRoughFwhm)
        self.log.logdebug("nSmallestSampleRoughFwhm: %d" % self.config.nSmallSampleRoughFwhm)

        # extracting the given number of most compact sources

        if True:
            # first pick n compact sources, then pick n bright sources of them
            if self.config.doUndistort:
                indicesSourcesSmall = numpy.argsort(fwhmUndistListForRoughFwhm)[:self.config.nSmallSampleRoughFwhm]
            else:
                indicesSourcesSmall = numpy.argsort(fwhmListForRoughFwhm)[:self.config.nSmallSampleRoughFwhm]
            magListForRoughFwhmSmall = magListForRoughFwhm[indicesSourcesSmall]
            indicesSourcesBrightOfSmall = numpy.argsort(magListForRoughFwhmSmall)[:self.config.nBrightSampleRoughFwhm]
            indicesSourcesPsfLike = indicesSourcesSmall[indicesSourcesBrightOfSmall]
        elif False:
            # first pick n brightest sources, then pick n compact sources of them
            indicesSourcesBright = numpy.argsort(magListForRoughFwhm)[:self.config.nBrightSampleRoughFwhm]
            if self.config.doUndistort:
                fwhmUndistListForRoughFwhmBright = fwhmUndistListForRoughFwhm[indicesSourcesBright]
                indicesSourcesSmallOfBright = numpy.argsort(fwhmUndistListForRoughFwhmBright)[:self.config.nSmallSampleRoughFwhm]
            else:
                fwhmListForRoughFwhmBright = fwhmListForRoughFwhm[indicesSourcesBright]
                indicesSourcesSmallOfBright = numpy.argsort(fwhmListForRoughFwhmBright)[:self.config.nSmallSampleRoughFwhm]
            indicesSourcesPsfLike = indicesSourcesBright[indicesSourcesSmallOfBright]
        else: # old; not optimized for edge CCDs where psf elongated and so psf sequence is broad in fwhm axis
            if self.config.doUndistort:
                indicesSourcesPsfLike = numpy.argsort(fwhmUndistListForRoughFwhm)[:self.config.nSmallSampleRoughFwhm]
            else:
                indicesSourcesPsfLike = numpy.argsort(fwhmListForRoughFwhm)[:self.config.nSmallSampleRoughFwhm]

        magListPsfLike = magListForRoughFwhm[indicesSourcesPsfLike]
        fwhmListPsfLike = fwhmListForRoughFwhm[indicesSourcesPsfLike]
        fwhmUndistListPsfLike = []
        fwhmRough = None
        fwhmUndistRough = None

        if self.config.statAlgRoughFwhm == "MEANCLIP":
            clipSigma = 3.0
            nIter = 3
            fwhmStat = afwMath.MEANCLIP
            sctrl = afwMath.StatisticsControl(clipSigma, nIter)
            stats = afwMath.makeStatistics(fwhmListPsfLike, fwhmStat, sctrl)
            fwhmRough = stats.getValue(fwhmStat)
            if self.config.doUndistort:
                stats = afwMath.makeStatistics(fwhmUndistListPsfLike, fwhmStat, sctrl)
                fwhmUndistRough = stats.getValue(fwhmStat)
        else:
            fwhmRough = numpy.median(fwhmListPsfLike)
            if self.config.doUndistort:
                fwhmUndistRough = numpy.median(fwhmUndistListPsfLike)

        if fwhmRough is None or numpy.isnan(fwhmRough) or fwhmRough <= 0:
            fwhmRough = -9999.0
        if fwhmUndistRough is None or numpy.isnan(fwhmUndistRough) or fwhmUndistRough <= 0:
            fwhmUndistRough = -9999.0

        data.magListPsfLike = magListPsfLike
        data.fwhmListPsfLike = fwhmListPsfLike
        data.fwhmRough = fwhmRough
        data.fwhmUndistListPsfLike = fwhmUndistListPsfLike
        data.fwhmUndistRough = fwhmUndistRough

        if self.debugFlag:
            print '*** fwhmListForRoughFwhm:', fwhmListForRoughFwhm
            print '*** fwhmRough:', fwhmRough
        self.log.info("fwhmRough: %f" % fwhmRough)
        self.metadata.set("fwhmRough", fwhmRough)

        if self.config.doPlots and dataRef is not None:
            #fig = plt.figure()
            fig = figure.Figure()
            canvas = FigCanvas(fig)
            pltMagFwhm = fig.add_subplot(1,1,1)
            pltMagFwhm.set_xlim(-20,-5)
            pltMagFwhm.set_ylim(0,20)
            pltMagFwhm.plot(data.magListAll, data.fwhmListAll, 'c+', label='all sample')
            pltMagFwhm.plot(magListForRoughFwhm, fwhmListForRoughFwhm, 'bx', label='sample for rough FWHM')
            pltMagFwhm.plot(magListPsfLike, fwhmListPsfLike, 'ko', label='coarse PSF-like sample')
            xx = [-20,-5]
            yy = [fwhmRough, fwhmRough]
            pltMagFwhm.plot(xx, yy, linestyle='dashed', label='tentative FWHM')
            pltMagFwhm.set_title('FWHM vs magnitudes', position=(0.5,1.05))
            pltMagFwhm.set_xlabel('magnitude (bottom:isnt / top:calib)')
            pltMagFwhm.set_ylabel('FWHM (pix)')
            pltMagFwhm.legend()
            try:
                zpFrame = exposure.getCalib().getMagnitude(1.0) # zp(mag/ADU/exptime)
            except: # photocal yet to be done or failure in photocal
                zpFrame = None
            if zpFrame and not numpy.isnan(zpFrame):
                pltMagFwhm2 = pltMagFwhm.twiny() # another x axis with the calibrated mag
                xlimCalib = [ x + zpFrame for x in pltMagFwhm.get_xlim() ]
                pltMagFwhm2.set_xlim(xlimCalib)
                pltMagFwhm2.set_ylim(pltMagFwhm.get_ylim())
                pltMagFwhm2.plot([-100,100],[-1,-1])

            fname = getFilename(dataRef, "plotSeeingRough")
            fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait',
                        papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

            del fig
            del pltMagFwhm
            del canvas
            if zpFrame and not numpy.isnan(zpFrame):
                del pltMagFwhm2
 
        return fwhmRough

    def getFwhmRobust(self, dataRef, data, fwhmRough, exposure):
        """Good Estimation of Final PSF"""

        # deriving final representative values of seeing, ellipticity and pa of elongation

        if data:
            fwhmMin = fwhmRough - self.config.fwhmMarginFinal
            fwhmMax = fwhmRough + self.config.fwhmMarginFinal
            nbin = (fwhmMax-fwhmMin) / self.config.fwhmBinSize
            # finding the mode
            histFwhm = numpy.histogram(data.fwhmListForPsfSeq, range=(fwhmMin, fwhmMax), bins=nbin)
            minval=0
            for i, num in enumerate(histFwhm[0]):
                if num > minval:
                    icand = i
                    minval = num

            numFwhmRobust = histFwhm[0][icand]
            fwhmRobust = histFwhm[1][icand] + 0.5*self.config.fwhmBinSize
            ellRobust = numpy.median(data.ellListPsfLikeRobust)
            ellPaRobust = numpy.median(data.ellPaListPsfLikeRobust)

            if fwhmRobust is None or numpy.isnan(fwhmRobust) or fwhmRobust <= 0:
                fwhmRobust = -9999.0
            if ellRobust is None or numpy.isnan(ellRobust) or ellRobust <= 0:
                ellRobust = -9999.0
            if ellPaRobust is None or numpy.isnan(ellPaRobust) or ellPaRobust < -90 or ellPaRobust > 90:
                ellPaRobust = -9999.0
        else:
            fwhmRobust = -9999.0
            ellRobust = -9999.0
            ellPaRobust = -9999.0

        self.log.info("Robust quantities: %f %f %f" % (fwhmRobust, ellRobust, ellPaRobust))
        self.metadata.set("fwhmRobust", fwhmRobust)
        self.metadata.set("ellRobust", ellRobust)
        self.metadata.set("ellPaRobust", ellPaRobust)

        if self.config.doPlots and dataRef is not None and data is not None:
            fig = figure.Figure()
            canvas = FigCanvas(fig)
            #fig = plt.figure()
            pltMagFwhm = fig.add_subplot(1,2,1)
            pltMagFwhm.set_xlim(-20,-5)
            pltMagFwhm.set_ylim(0,20)
            pltMagFwhm.plot(data.magListAll, data.fwhmListAll, 'c+', label='all sample')
            #pltMagFwhm.plot(data.magListPsfLike, data.fwhmListPsfLike, 'bx', label='coarse PSF-like sample')
            pltMagFwhm.plot(data.magListForRoughFwhm, data.fwhmListForRoughFwhm, 'bx', label='sample for rough FWHM')
            #pltMagFwhm.plot(data.magListPsfLike, data.fwhmListPsfLike, 'bx', label='sample for rough fwhm')
            pltMagFwhm.plot(data.magListPsfLikeRobust, data.fwhmListPsfLikeRobust, 'mo',
                            label='best-effort PSF-like sample')
            xx = [-20,-5]; yy = [fwhmRobust, fwhmRobust]
            yy_median = [data.medianFwhmPsfSeq, data.medianFwhmPsfSeq]
            yy_sigma_h = [data.fwhmPsfSeqMax, data.fwhmPsfSeqMax]
            yy_sigma_l = [data.fwhmPsfSeqMin, data.fwhmPsfSeqMin]
            pltMagFwhm.plot(xx, yy, linestyle='dashed', color='blue', label='resultant FWHM: %4.2f' % fwhmRobust)
            pltMagFwhm.plot(xx, yy_median, linestyle='dashed', color='red', label='median FHWM of sequence')
            pltMagFwhm.plot(xx, yy_sigma_l, linestyle='-.', color='blue', label='sigma range')
            pltMagFwhm.plot(xx, yy_sigma_h, linestyle='-.', color='blue', label=None)
            pltMagFwhm.set_title('FWHM vs magnitudes', position=(0.45,1.05))
            pltMagFwhm.set_xlabel('magnitude (bottom:isnt / top:calib)')
            pltMagFwhm.set_ylabel('FWHM (pix)')
            pltMagFwhm.legend()
            try:
                zpFrame = exposure.getCalib().getMagnitude(1.0) # zp(mag/ADU/exptime)
            except: # photocal yet to be done or failure in photocal
                zpFrame = None
            if zpFrame and not numpy.isnan(zpFrame):
                pltMagFwhm2 = pltMagFwhm.twiny() # another x axis with the calibrated mag
                xlimCalib = [ x + zpFrame for x in pltMagFwhm.get_xlim() ]
                pltMagFwhm2.set_xlim(xlimCalib)
                pltMagFwhm2.set_ylim(pltMagFwhm.get_ylim())
                pltMagFwhm2.plot([-100,100],[-1,-1])

            pltHistFwhm = fig.add_subplot(1,2,2)
            pltHistFwhm.hist(data.fwhmListForPsfSeq, range=(fwhmMin, fwhmMax), bins=nbin, orientation='horizontal', histtype='bar')
            pltHistFwhm.plot([0, numFwhmRobust*1.2], [fwhmRobust, fwhmRobust], linestyle='dashed', color='black', label=None)
            pltHistFwhm.set_title('histogram of FWHM')
            pltHistFwhm.set_ylabel('FWHM of best-effort PSF-like sources')
            pltHistFwhm.set_xlabel('number of sources')
            pltHistFwhm.legend()

            fname = getFilename(dataRef, "plotSeeingRobust")
            fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait',
                        papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

            del fig
            del pltMagFwhm
            del pltHistFwhm
            del canvas
            if zpFrame and not numpy.isnan(zpFrame):
                del pltMagFwhm2

        data.fwhmRobust = fwhmRobust
        data.ellRobust = ellRobust
        data.ellPaRobust = ellPaRobust

        return data

    def getStarCandidateList(self, dataRef, data, fwhmRough, magLimPsfSeq):
        """Good Estimation of Final PSF sample"""

        # deriving sigma of scatter of the psf-like sequence around the fwhmRough in the 'mag vs fwhm' space
        # sample around the fwhmRough is limited by the given config fwhmMarginFinal [pix]

        fwhmMin = fwhmRough - self.config.fwhmMarginFinal
        fwhmMax = fwhmRough + self.config.fwhmMarginFinal

        if self.config.doUndistort:
            fwhmListAllToEval = data.fwhmUndistListAll
        else:
            fwhmListAllToEval = data.fwhmListAll
        indicesSourcesForPsfSeq = [ i for i in data.indicesSourcesFwhmRange if
                                    fwhmListAllToEval[i] > fwhmMin and fwhmListAllToEval[i] < fwhmMax and data.magListAll[i] < magLimPsfSeq ]
        magListForPsfSeq = data.magListAll[indicesSourcesForPsfSeq]
        fwhmListForPsfSeq = data.fwhmListAll[indicesSourcesForPsfSeq]
        if self.config.doUndistort:
            fwhmUndistListForPsfSeq = data.fwhmUndistListAll[indicesSourcesForPsfSeq]
            fwhmListForPsfSeqToEval = fwhmUndistListForPsfSeq
        else:
            fwhmUndistListForPsfSeq = []
            fwhmListForPsfSeqToEval = fwhmListForPsfSeq

        #data.indicesForPsfSeq = indicesSourcesForPsfSeq
        #data.magListForPsfSeq = magListForPsfSeq
        #data.fwhmListForPsfSeq = fwhmListForPsfSeq

        # deriving fwhwRobust and sigma
        clipSigma = self.config.psfSeqStatNsigma
        nIter = self.config.psfSeqStatNiter
        fwhmStat = afwMath.MEDIAN
        sigmaStat = afwMath.STDEVCLIP
        sctrl = afwMath.StatisticsControl(clipSigma, nIter)

        if len(fwhmListForPsfSeqToEval) > 0:
            stats = afwMath.makeStatistics(fwhmListForPsfSeqToEval, fwhmStat | sigmaStat, sctrl)
            medianFwhmPsfSeq = stats.getValue(fwhmStat)
            sigmaFwhmPsfSeq = stats.getValue(sigmaStat)

            # extracting psf candidates by limting given number of sigma around the psf sequence
            # maglim is extended toward the faint end if desired
            magLimPsfCand = magLimPsfSeq + self.config.magLimitFaintExtension

            fwhmPsfSeqMin = medianFwhmPsfSeq - self.config.fwhmMarginNsigma*sigmaFwhmPsfSeq
            fwhmPsfSeqMax = medianFwhmPsfSeq + self.config.fwhmMarginNsigma*sigmaFwhmPsfSeq

            indicesSourcesPsfLikeRobust = [ i for i in data.indicesSourcesFwhmRange if
                                            fwhmListAllToEval[i] > fwhmPsfSeqMin and fwhmListAllToEval[i] < fwhmPsfSeqMax 
                                            and data.magListAll[i] < magLimPsfCand  ]

            numFwhmPsfLikeRobust = len(indicesSourcesPsfLikeRobust)
            if numFwhmPsfLikeRobust < 1:
                self.log.warn("No sources selected in robust seeing estimation")
                return None
        else:
            self.log.warn("No sources selected in psf sequence candidates for robust seeing estimation")
            return None

        print '*** medianFwhmPsfSeq:', medianFwhmPsfSeq
        self.metadata.set("medianFwhmPsfSeq", medianFwhmPsfSeq)
        print '*** sigmaFwhmPsfSeq:', sigmaFwhmPsfSeq
        self.metadata.set("sigmaFwhmPsfSeq", sigmaFwhmPsfSeq)
        print '*** fwhmPsfSequence(min, max): (%5.3f %5.3f) (pix)' % (fwhmPsfSeqMin, fwhmPsfSeqMax)
        self.metadata.set("minFwhmPsfSeq", fwhmPsfSeqMin)
        self.metadata.set("maxFwhmPsfSeq", fwhmPsfSeqMax)
        print '*** numFwhmPsfLikeRobust:', numFwhmPsfLikeRobust
        self.metadata.set("numFwhmPsfLikeRobust", numFwhmPsfLikeRobust)

        magListPsfLikeRobust = data.magListAll[indicesSourcesPsfLikeRobust]
        fwhmListPsfLikeRobust = data.fwhmListAll[indicesSourcesPsfLikeRobust]
        ellListPsfLikeRobust = data.ellListAll[indicesSourcesPsfLikeRobust]
        ellPaListPsfLikeRobust = data.ellPaListAll[indicesSourcesPsfLikeRobust]
        AEllListPsfLikeRobust = data.AEllListAll[indicesSourcesPsfLikeRobust]
        BEllListPsfLikeRobust = data.BEllListAll[indicesSourcesPsfLikeRobust]

        elle1e2ListPsfLikeRobust = data.elle1e2ListAll[indicesSourcesPsfLikeRobust]
        e1ListPsfLikeRobust = data.e1ListAll[indicesSourcesPsfLikeRobust]
        e2ListPsfLikeRobust = data.e2ListAll[indicesSourcesPsfLikeRobust]

        xListPsfLikeRobust = data.xListAll[indicesSourcesPsfLikeRobust]
        yListPsfLikeRobust = data.yListAll[indicesSourcesPsfLikeRobust]
        IxxListPsfLikeRobust = data.IxxListAll[indicesSourcesPsfLikeRobust]
        IyyListPsfLikeRobust = data.IyyListAll[indicesSourcesPsfLikeRobust]
        IxyListPsfLikeRobust = data.IxyListAll[indicesSourcesPsfLikeRobust]

        data.indicesSourcesPsfLikeRobust = indicesSourcesPsfLikeRobust
        data.magListPsfLikeRobust = magListPsfLikeRobust
        data.fwhmListPsfLikeRobust = fwhmListPsfLikeRobust
        data.ellListPsfLikeRobust = ellListPsfLikeRobust
        data.ellPaListPsfLikeRobust = ellPaListPsfLikeRobust
        data.AEllListPsfLikeRobust = AEllListPsfLikeRobust
        data.BEllListPsfLikeRobust = BEllListPsfLikeRobust

        data.elle1e2ListPsfLikeRobust = elle1e2ListPsfLikeRobust
        data.e1ListPsfLikeRobust = e1ListPsfLikeRobust
        data.e2ListPsfLikeRobust = e2ListPsfLikeRobust

        data.xListPsfLikeRobust = xListPsfLikeRobust
        data.yListPsfLikeRobust = yListPsfLikeRobust
        data.IxxListPsfLikeRobust = IxxListPsfLikeRobust
        data.IyyListPsfLikeRobust = IyyListPsfLikeRobust
        data.IxyListPsfLikeRobust = IxyListPsfLikeRobust
        data.indicesForPsfSeq = indicesSourcesForPsfSeq
        data.magListForPsfSeq = magListForPsfSeq
        data.fwhmListForPsfSeq = fwhmListForPsfSeq
        data.fwhmUndistListForPsfSeq = fwhmUndistListForPsfSeq
        data.magLimPsfSeq = magLimPsfSeq
        data.magLimPsfCand = magLimPsfCand
        data.medianFwhmPsfSeq = medianFwhmPsfSeq
        data.sigmaFwhmPsfSeq = sigmaFwhmPsfSeq
        data.fwhmPsfSeqMin = fwhmPsfSeqMin
        data.fwhmPsfSeqMax = fwhmPsfSeqMax

        return data

#measAlg.starSelectorRegistry.register("mitaka", SizeMagnitudeMitakaStarSelector)

def getFilename(dataRef, dataset):
    fname = dataRef.get(dataset + "_filename")[0]
    directory = os.path.dirname(fname)
    if not os.path.exists(directory):
        try: # handling a possible race condition
            os.makedirs(directory)
        except OSError, e:
            if e.errno != errno.EEXIST:
                raise
        if not os.path.exists(directory):
            raise RuntimeError, "Unable to create output directory for qa output directory '%s'" % (directory)

    return fname


measAlg.starSelectorRegistry.register("mitaka", SizeMagnitudeMitakaStarSelector)
