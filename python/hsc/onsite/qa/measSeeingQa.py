#!/usr/bin/env python

"""
TODO:
 This module is to use the same algorithms to select stars & measure seeing
  as in the SC-RCM running for Suprime-Cam observation (publised in
  Furusawa+2011,PASJ,63,S581).
  Implemented for onsite QA output
"""

import sys
import os, re
import math
import lsst.pex.config as pexConfig
import lsst.pex.logging as pexLog
#import lsst.afw.detection as afwDetection
#import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.table as afwTable

import numpy

# this order is necessary to avoid X connection from this script
import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
from matplotlib import patches as patches

##plt.switch_backend('Agg') # experimental option to work around unknown conflict in the Agg backend 
### Be sure that other python modules do not call matplotlib.pyplot with TkAgg or non-Agg backend before this module.

class QaSeeingConfig(pexConfig.Config):
    fwhmIni = pexConfig.Field(
        dtype = float,
        doc = 'Inintial guess of seeing (pix)',
        default = 3.465,
        )
    fracSrcIni = pexConfig.Field(
        dtype = float, 
        doc = 'What fraction of sources from the brightest is to be included for initial guess of seeing to avoid cosmic rays which dominate faint magnitudes', 
        default = 0.15, # is good for SC with 5-sigma detection
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
    nSampleRoughFwhm = pexConfig.Field(
        dtype = int,
        doc = 'Number of smallest objects which are used to determine rough-interim seeing',
        default = 30,
        )
    fwhmMarginFinal = pexConfig.Field(
        dtype = float,
        doc = 'How many pixels around the peak are used for final seeing estimation as mode',
        default = 1.5,
        ) 
    camera = pexConfig.Field(
        dtype = str, 
        doc = 'camera same as in root.camera', 
        default = 'hsc',
        )
   
def measureSeeingQaTest(exposure, config):
    print '**** measSeeing: qa.seeing.fwhmIni: ', config.qa.seeing.fwhmIni
    print '**** measSeeing: qa.seeing.fwhmMin: ', config.qa.seeing.fwhmMin
    print '**** measSeeing: qa.seeing.fwhmMax: ', config.qa.seeing.fwhmMax
    print '**** measSeeing: qa.flatness.meshX: ', config.qa.flatness.meshX
    print '**** config all:'
    print config
    
def measureSeeingQa(exposure, catalog, config, dataId, debugFlag=False, plotFlag=True, plotbasedir=None, butler=None, log=None):
    """
    This function measures seeing (mode) and ellipticity (median) representative of the input frame
    by using the input source set.

    Requirements:
       This function has optimized for reasonably deep detection. Please use 2 sigma detection for
       the input sourceset. 

       exposure is used for metadata filling
       
    """
    if log is None:
        log = pexLog.Log.getDefaultLog()

    if config is not None:
        seeingConfig = config.qa.seeing
        fwhmIni = seeingConfig.fwhmIni  # initial guess of seeing size (pix)
        fwhmMin = seeingConfig.fwhmMin  # minimum seeing to treat (pix)
        fwhmMax = seeingConfig.fwhmMax  # maximum seeing to treat (pix)
        nbinMagHist = seeingConfig.nbinMagHist
        magMinHist = seeingConfig.magMinHist
        magMaxHist = seeingConfig.magMaxHist
        nSampleRoughFwhm = seeingConfig.nSampleRoughFwhm
        fwhmMarginFinal = seeingConfig.fwhmMarginFinal
        fracSrcIni = seeingConfig.fracSrcIni
        
    else:
        fwhmIni = 0.7/0.202  # pix <- assuming SC. 
        fwhmMin = 1.5        # pix
        fwhmMax = 12.0       # pix

        #nbinMagHist = 10
        nbinMagHist= 20 / 0.25
        magMinHist = -20.0
        magMaxHist = 0.0

        #nSampleRoughFwhm = 50 # default objects
        nSampleRoughFwhm = 30 # objects
        fwhmMarginFinal = 1.5  # pix

        fracSrcIni = 0.15

    # -- Filtering sources with rough min-max fwhm
    xListAll = []
    yListAll = []
    magListAll = []
    fwhmListAll = []
    ellListAll = []
    ellPaListAll = []
    IxxListAll = []
    IyyListAll = []
    IxyListAll = []
    AEllListAll = []
    BEllListAll = []        
    objFlagListAll = []
    indicesSourcesFwhmRange = [] # indices of sources in acceptable fwhm range 

#    flags = maUtils.getDetectionFlags()
#    psfLikeFlag = flags['STAR']
    catalogSchema = catalog.table.getSchema()
    nameSaturationFlag = 'flags.pixel.saturated.any'
    nameSaturationCenterFlag = 'flags.pixel.saturated.center'
    keySaturationFlag = catalogSchema.find(nameSaturationFlag).key
    keySaturationCenterFlag = catalogSchema.find(nameSaturationCenterFlag).key    

    print 'keySaturationFlag:',  keySaturationFlag
    print 'keySaturationCenterFlag:',  keySaturationCenterFlag    
    iGoodId = 0
    for iseq, source in enumerate(catalog):

        xc = source.getX()   
        yc = source.getY()   
        Ixx = source.getIxx() # luminosity-weighted 2nd moment of pixels
        Iyy = source.getIyy() 
        Ixy = source.getIxy()
        sigma = math.sqrt(0.5*(Ixx+Iyy))
        fwhm = 2.*math.sqrt(2.*math.log(2.)) * sigma # (pix) assuming Gaussian

        saturationFlag = source.get(keySaturationFlag)
        saturationCenterFlag = source.get(keySaturationCenterFlag)        
#        objFlag = source.getFlagForDetection()
        isSaturated = (saturationFlag | saturationCenterFlag)        
        if debugFlag:
            print 'objFlag SatAny %x  SatCenter %x  isSaturated %x' % (saturationFlag, saturationCenterFlag, isSaturated)

        fluxAper = source.getApFlux()
        fluxErrAper =  source.getApFluxErr() 
        fluxGauss = source.getModelFlux()
        fluxErrGauss =  source.getModelFluxErr()
        fluxPsf = source.getPsfFlux()
        fluxErrPsf = source.getPsfFluxErr()

        # which flux do you use:
        fluxForSeeing = fluxAper 

        # validity check
#        if math.isnan(fwhm) or math.isnan(Ixx) or math.isnan(fluxForSeeing) or fwhm <= 0 or fluxForSeeing <= 0 or (Ixx == 0 and Iyy == 0 and Ixy == 0) or (objFlag & saturationFlag) > 0:
        if math.isnan(fwhm) or math.isnan(Ixx) or math.isnan(fluxForSeeing) or fwhm <= 0 or fluxForSeeing <= 0 or (Ixx == 0 and Iyy == 0 and Ixy == 0) or isSaturated is True:
            # this sample is excluded
            continue

        mag = -2.5*numpy.log10(fluxForSeeing)
        magListAll.append(mag)
        fwhmListAll.append(fwhm)
        xListAll.append(xc)
        yListAll.append(yc)
        IxxListAll.append(Ixx)
        IyyListAll.append(Iyy)
        IxyListAll.append(Ixy)

        if True: # definition by SExtractor
            Val1 = 0.5*(Ixx+Iyy)
            Ixx_Iyy = Ixx-Iyy
            Val2 = 0.25*Ixx_Iyy*Ixx_Iyy + Ixy*Ixy
            if Val2 >= 0 and (Val1-math.sqrt(Val2)) > 0:
                aa = math.sqrt( Val1 + math.sqrt(Val2) )
                bb = math.sqrt( Val1 - math.sqrt(Val2) )
                ell =  1. - bb/aa
                if math.fabs(Ixx_Iyy) > 1.0e-10:
                    ellPa = 0.5 * math.degrees(math.atan(2*Ixy / math.fabs(Ixx_Iyy)))
                else:
                    ellPa = 0.0
            else:
                ell = None
                ellPa = None
                aa = None
                bb = None
        else: # definition by Kaiser
            # e=sqrt(e1^2+e2^2) where e1=(Ixx-Iyy)/(Ixx+Iyy), e2=2Ixy/(Ixx+Iy)
            # SExtractor's B/A=sqrt((1-e)/(1+e)), ell=1-B/A
            e1 = (Ixx-Iyy)/(Ixx+Iyy)
            if e1 > 0: 
                e2 = 2.0*Ixy/(Ixx+Iyy)
                ell = math.sqrt(e1*e1 + e2*e2)
                fabs_Ixx_Iyy = math.fabs(Ixx-Iyy)
                if fabs_Ixx_Iyy > 1.0e-10:
                    ellPa = 0.5 * math.degrees(math.atan(2*Ixy / fabs_Ixx_Iyy))
                else:
                    ellPa = 0.0
            else:
                ell = None
                ellPa = None

        if ellPa is not None:
            ellPa = 90. - ellPa ## definition of PA to be confirmed

        if debugFlag:
            print '*** %d : Ixx: %f Iyy: %f Ixy: %f fwhm: %f flux: %f isSatur: %x' % (iseq, Ixx, Iyy, Ixy, fwhm, fluxAper, isSaturated)

        ellListAll.append( ell )
        ellPaListAll.append( ellPa )
        AEllListAll.append( aa ) # sigma in long axis
        BEllListAll.append( bb ) # sigma in short axis
        #objFlagListAll.append(objFlag)
        
        if fwhm > fwhmMin and fwhm < fwhmMax:
            indicesSourcesFwhmRange.append(iGoodId)

        # this sample is added to the good sample
        iGoodId += 1

    #import pdb; pdb.set_trace()

    # -- Deriving the normalized cumulative magnitude histogram
    magListAll = numpy.array(magListAll)
    fwhmListAll = numpy.array(fwhmListAll)
    ellListAll = numpy.array(ellListAll)
    ellPaListAll = numpy.array(ellPaListAll)    
    AEllListAll = numpy.array(AEllListAll)
    BEllListAll = numpy.array(BEllListAll)    
    xListAll = numpy.array(xListAll)
    yListAll = numpy.array(yListAll)    
    IxxListAll = numpy.array(IxxListAll)    
    IyyListAll = numpy.array(IyyListAll)    
    IxyListAll = numpy.array(IxyListAll)    

    #indicesSourcesFwhmRange = numpy.array(indicesSourcesFwhmRange)
    magListFwhmRange = magListAll[indicesSourcesFwhmRange] # magList only for sources within the FWHM range
    magHist = numpy.histogram(magListFwhmRange, range=(magMinHist,magMaxHist), bins=nbinMagHist)  # n.b. magHist[1] is sorted so index=0 is brightest

    sumAll = magHist[0].sum()
    #print '*** sumAll: ', sumAll
    log.log(log.INFO, "QaSeeing: total number of objects in the first dataset (sumAll): %d" % sumAll)

    magCumHist = list(magHist)
    magCumHist[0] = (magCumHist[0].cumsum())
    #print '------ mag cum hist ------'
    log.log(log.INFO, "QaSeeing: calculating cumulative n(m)")    
#    print magCumHist
    magCumHist[0] = magCumHist[0] / float(sumAll)
    #print '------ normalized mag cum hist ------'
    log.log(log.INFO, "QaSeeing: normalizing cumulative n(m) by the sumAll")        
#    print magCumHist

    # -- Estimating mag limit based no the cumlative mag histogram

    magLim = None
    for i, cumFraction in enumerate(magCumHist[0]):
        if cumFraction >= fracSrcIni:
            magLim = magCumHist[1][i] # magLim is the mag which exceeds the cumulative n(m) of 0.15
            break
    if not magLim:
        log.log(log.WARN, "QaSeeing: Error: cumulative magnitude histogram does not exceed 0.15.")
        return None

    #print '*** magLim: ', magLim
    log.log(log.INFO, "QaSeeing: mag limit auto-dertermined: %5.2f or %f (ADU)" % (magLim, numpy.power(10, -0.4*magLim)))

    if debugFlag:
        log.log(log.DEBUG, "QaSeeing: magHist: %s" % magHist)
        log.log(log.DEBUG, "QaSeeing: cummurative magHist: %s" % magCumHist)
        #print magHist
        #print magCumHist

    # making seeing-debug1 plots
    if plotFlag is True: 
        log.log(log.INFO, "QaSeeing: Entered plotting")
        frameId = exposure.getMetadata().get('FRAMEID')

        qaOutputDirName = os.path.dirname(butler.get('src_filename', dataId)[0])
        #if os.path.exists(plotbasedir):
        #    qaOutputDirName = plotbasedir
        #elif butler is not None:
        #    qaOutputDirName = os.path.dirname(butler.get('src_filename', dataId)[0])
        #else:
        #    qaOutputDirName = os.path.abspath('.')
        #log.log(log.INFO, "QaSeeing: Qa plots are stored in: %s", qaOutputDirName)            
            
        # making a figure
        fig = plt.figure()
        pltMagHist = fig.add_subplot(2,1,1)
        pltMagHist.hist(magListFwhmRange, bins=magHist[1], orientation='vertical')
        pltMagHist.set_title('histogram of magnitudes')
        #        pltMagHist.set_xlabel('magnitude instrumental')
        pltMagHist.set_ylabel('number of samples')
        pltMagHist.legend()

        pltCumHist = fig.add_subplot(2,1,2)
        # histtype=bar,barstacked,step,stepfilled
        pltCumHist.hist(magListFwhmRange, bins=magCumHist[1], normed=True, cumulative=True, orientation='vertical', histtype='step')
        xx = [magLim, magLim]; yy = [0, 1]
        pltCumHist.plot(xx, yy, linestyle='dashed', label='mag limit') # solid,dashed,dashdot,dotted        
        pltCumHist.set_title('cumulative histogram of magnitudes')
        pltCumHist.set_xlabel('magnitude instrumental')
        pltCumHist.set_ylabel('Nsample scaled to unity')
        pltCumHist.legend()
        #        plt.draw()
        #        plt.show()

        fname = os.path.join(qaOutputDirName, 'seeing1_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltMagHist
        del pltCumHist        

    # -- Estimating roughly-estimated FWHM for sources with mag < magLim

    # Do we use STAR-Like flag in detection?
    useStarFlag = False
    if useStarFlag:
        #flags = maUtils.getDetectionFlags()
        #psfLikeFlag = flags['STAR']
        #psfLikeFlag = flags['BINNED1'] | flags['STAR']
        indicesSourcesPsfLike = \
                              [i for i in indicesSourcesFwhmRange if ((magListAll[i]<magLim) and (objFlagListAll[i] & psfLikeFlag))]
    else:
        indicesSourcesPsfLike = [i for i in indicesSourcesFwhmRange if magListAll[i]<magLim ]        

    magListPsfLike = magListAll[indicesSourcesPsfLike]
    fwhmListPsfLike = fwhmListAll[indicesSourcesPsfLike]
    ellListPsfLike = ellListAll[indicesSourcesPsfLike]
    ellPaListPsfLike = ellPaListAll[indicesSourcesPsfLike]
    AEllListPsfLike = AEllListAll[indicesSourcesPsfLike]
    BEllListPsfLike = BEllListAll[indicesSourcesPsfLike]
    xListPsfLike = xListAll[indicesSourcesPsfLike]
    yListPsfLike = yListAll[indicesSourcesPsfLike]
    IxxListPsfLike = IxxListAll[indicesSourcesPsfLike]
    IyyListPsfLike = IyyListAll[indicesSourcesPsfLike]
    IxyListPsfLike = IxyListAll[indicesSourcesPsfLike]
    
    #print '*** indicesSourcesFwhmRange:', indicesSourcesFwhmRange
    #print '*** indicesSoucresPsfLike:',  indicesSourcesPsfLike
    #print '*** magListPsfLike:',  magListPsfLike
    #print '*** fwhmListPsfLike:',  fwhmListPsfLike

    log.log(log.INFO, "QaSeeing: nSampleRoughFwhm: %d" % nSampleRoughFwhm)

    fwhmListForRoughFwhm = numpy.sort(fwhmListPsfLike)[:nSampleRoughFwhm] # fwhmListPsfLike is not changed with numpy.sort()
    fwhmRough = numpy.median(fwhmListForRoughFwhm)

    if debugFlag:
        print '*** fwhmListPsfLike:', fwhmListPsfLike
        print '*** fwhmListForRoughFwhm:', fwhmListForRoughFwhm
        print '*** fwhmRough:', fwhmRough
    log.log(log.INFO, "QaSeeing: fwhmRough: %f" % fwhmRough)
    
    # making seeing-debug2 plots
    if plotFlag: 
        fig = plt.figure()
        pltMagFwhm = fig.add_subplot(1,1,1)
        pltMagFwhm.set_xlim(-20,-5)
        pltMagFwhm.set_ylim(0,20)
        pltMagFwhm.plot(magListPsfLike, fwhmListPsfLike, 'bx', label='coarse PSF-like sample')
        xx = [-20,-5]; yy = [fwhmRough, fwhmRough]
        pltMagFwhm.plot(xx, yy, linestyle='dashed', label='tentative FWHM') # solid,dashed,dashdot,dotted
        pltMagFwhm.set_title('FWHM vs magnitudes')
        pltMagFwhm.set_xlabel('magnitude instrumental')
        pltMagFwhm.set_ylabel('FWHM (pix)')

        pltMagFwhm.legend()
#        plt.show()
#        plt.draw()
        fname = os.path.join(qaOutputDirName, 'seeing2_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltMagFwhm


    # ---- From Here, Good Estimation of Final PSF

    fwhmMin = fwhmRough - fwhmMarginFinal
    fwhmMax = fwhmRough + fwhmMarginFinal
##    fwhmMin = fwhmRough - 1.0 # pix
##    fwhmMax = fwhmRough + 1.0 # pix 
##    fwhmMin = fwhmRough - 0.5 # pix
##    fwhmMax = fwhmRough + 0.5 # pix 

    # Do we use STAR-Like flag in detection? Do we use Magnitude Limit?
    useStarFlag = False
    if useStarFlag:
        #flags = maUtils.getDetectionFlags()
        psfLikeFlag = flags['STAR']
        #psfLikeFlag = flags['BINNED1'] | flags['STAR']
        indicesSourcesPsfLikeRobust = [ \
            i for i in indicesSourcesFwhmRange \
            if ((fwhmListAll[i]>fwhmMin) and (fwhmListAll[i]<fwhmMax) \
                                              and (objFlagListAll[i] & psfLikeFlag)) \
            ]
    else:
        indicesSourcesPsfLikeRobust = [ \
            i for i in indicesSourcesFwhmRange \
            if ((fwhmListAll[i]>fwhmMin) and (fwhmListAll[i]<fwhmMax)) \
            ] 

    # error if too few sample remain
    if len(indicesSourcesPsfLikeRobust) < 1:
        #print '*** Number of sample stars are too small. None is returned.'
        log.log(log.INFO, "QaSeeing: number of sample stars are too small. None is returned.")
        magListPsfLikeRobust = None
        fwhmListPsfLikeRobust = None
        ellListPsfLikeRobust = None
        ellPaListPsfLikeRobust = None        
        AEllListPsfLikeRobust = None
        BEllListPsfLikeRobust = None        
        fwhmRobust = None
        ellRobust = None
        ellPaRobust = None
        xListPsfLikeRobust = None
        yListPsfLikeRobust =  None
        IxxListPsfLikeRobust =  None
        IyyListPsfLikeRobust =  None
        IxyListPsfLikeRobust =  None
    else:
        magListPsfLikeRobust = magListAll[indicesSourcesPsfLikeRobust]
        fwhmListPsfLikeRobust = fwhmListAll[indicesSourcesPsfLikeRobust]
        ellListPsfLikeRobust = ellListAll[indicesSourcesPsfLikeRobust]
        ellPaListPsfLikeRobust = ellPaListAll[indicesSourcesPsfLikeRobust]
        AEllListPsfLikeRobust = AEllListAll[indicesSourcesPsfLikeRobust]
        BEllListPsfLikeRobust = BEllListAll[indicesSourcesPsfLikeRobust]
        xListPsfLikeRobust = xListAll[indicesSourcesPsfLikeRobust]
        yListPsfLikeRobust = yListAll[indicesSourcesPsfLikeRobust]
        IxxListPsfLikeRobust = IxxListAll[indicesSourcesPsfLikeRobust]
        IyyListPsfLikeRobust = IyyListAll[indicesSourcesPsfLikeRobust]
        IxyListPsfLikeRobust = IxyListAll[indicesSourcesPsfLikeRobust]
        
        #print '*** fwhmListPsfLikeRobust:', fwhmListPsfLikeRobust
        log.log(log.INFO, "QaSeeing: fwhmListPsfLikeRobust: %s" % fwhmListPsfLikeRobust)
        
        binwidth = 0.2 
        nbin = (fwhmMax-fwhmMin)/binwidth
        # finding the mode
        histFwhm = numpy.histogram(fwhmListPsfLikeRobust, bins=nbin)
        minval=0
        for i, num in enumerate(histFwhm[0]):
            if num > minval:
                icand = i
                minval = num

        fwhmRobust = histFwhm[1][icand]
        ellRobust = numpy.median(ellListPsfLikeRobust)
        ellPaRobust = numpy.median(ellPaListPsfLikeRobust) 

    # making seeing-debug3 plots
    if plotFlag: 
        fig = plt.figure()
        pltMagFwhm = fig.add_subplot(1,2,1)
        pltMagFwhm.set_xlim(-20,-5)
        pltMagFwhm.set_ylim(0,20)
        pltMagFwhm.plot(magListAll, fwhmListAll, '+', label='all sample')
        pltMagFwhm.plot(magListPsfLike, fwhmListPsfLike, 'bx', label='coarse PSF-like sample')
        pltMagFwhm.plot(magListPsfLikeRobust, fwhmListPsfLikeRobust, 'mo', label='best-effort PSF-like sample')        
        xx = [-20,-5]; yy = [fwhmRobust, fwhmRobust]
        pltMagFwhm.plot(xx, yy, linestyle='dashed', label='resultant FWHM') # solid,dashed,dashdot,dotted
        pltMagFwhm.set_title('FWHM vs magnitudes')
        pltMagFwhm.set_xlabel('magnitude instrumental')
        pltMagFwhm.set_ylabel('FWHM (pix)')
#        plt.draw()
        pltMagFwhm.legend()

        pltHistFwhm = fig.add_subplot(1,2,2)
        pltHistFwhm.hist(fwhmListPsfLikeRobust, bins=nbin, orientation='horizontal', histtype='bar') 
        pltHistFwhm.set_title('histogram of FWHM')
        pltHistFwhm.set_ylabel('FWHM of best-effort PSF-like sources')
        pltHistFwhm.set_xlabel('number of sources')
        pltHistFwhm.legend()

#        plt.show()
        fname = os.path.join(qaOutputDirName, 'seeing3_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltHistFwhm

    #================================================================================================
    # making fwhm-map plots
    if True:
        xSize=2048.0
        ySize=4177.0
        facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
        wFig = xSize * facSize * 1.3
        hFig = ySize * facSize
        fig = plt.figure(figsize=(wFig,hFig))
        pltFwhmMap = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height
        
#        fig = plt.figure(1, figsize=(2048,4177))
#        pltFwhmMap = fig.add_subplot(1,1,1)
        pltFwhmMap.set_xlim(0, xSize)
        pltFwhmMap.set_ylim(0, ySize)
        pointSize = math.pi*(10*fwhmListPsfLikeRobust/2)**2. # 10pix=2arcsec fwhm = 50 point radius
#        pltFwhmMap.plot(xListPsfLikeRobust, yListPsfLikeRobust, markersize=pointSize, marker='o', label='PSF sample')
        pltFwhmMap.scatter(xListPsfLikeRobust, yListPsfLikeRobust, s=pointSize, marker='o', color=None, facecolor=(1,1,1,0), label='PSF sample')
        ###pltFwhmMap.legend()

        # reference sample point
        fwhmPix = numpy.array([5.0]) # pixel in fwhm
        pointSize = math.pi*(10*fwhmPix/2)**2.
        pltFwhmMap.scatter([0.1*xSize], [0.9*ySize], s=pointSize, marker='o', color='magenta', facecolor=(1,1,1,0), label='PSF sample')        
        plt.text(0.1 * xSize, 0.9 * ySize, 'fwhm=%4.1f pix' % fwhmPix, ha='center', va='top')

        pltFwhmMap.set_title('FWHM of PSF sources')
        pltFwhmMap.set_xlabel('X (pix)')
        pltFwhmMap.set_ylabel('Y (pix)')
#        plt.draw()
#        plt.show()
        fname = os.path.join(qaOutputDirName, 'fwhmmap_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltFwhmMap


    #================================================================================================
    # making psfellipse-map plots
    if True:
        xSize=2048.0
        ySize=4177.0
        facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
        wFig = xSize * facSize * 1.2
        hFig = ySize * facSize
        fig = plt.figure(figsize=(wFig,hFig))
        pltEllipseMap = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height
        
#        fig = plt.figure(1, figsize=(2048,4177))
#        pltEllipseMap = fig.add_subplot(1,1,1)
        pltEllipseMap.set_xlim(0, xSize)
        pltEllipseMap.set_ylim(0, ySize)
#        pltEllipseMap.plot(xListPsfLikeRobust, yListPsfLikeRobust, markersize=pointSize, marker='o', label='PSF sample')
        scaleFactor = 100.  # 10pix~2arcsec (width) -> 1000pix
        for xEllipse, yEllipse, aa, bb, ellPa in zip(xListPsfLikeRobust, yListPsfLikeRobust, AEllListPsfLikeRobust, BEllListPsfLikeRobust, ellPaListPsfLikeRobust):
            if all([aa, bb, ellPa]):
                ell = patches.Ellipse((xEllipse, yEllipse), 2.*aa*scaleFactor, 2.*bb*scaleFactor, angle=ellPa, linewidth=2., fill=False, zorder=2)
            pltEllipseMap.add_patch(ell)

        # reference sample point
        fwhmPix = aa = bb = 2.5 # pixel in A, B (in half width)
        ell = patches.Ellipse((0.1*xSize, 0.9*ySize), 2.*aa*scaleFactor, 2.*bb*scaleFactor, angle=0., linewidth=4., color='magenta', fill=False, zorder=2)
        pltEllipseMap.add_patch(ell)
        plt.text(0.1 * xSize, 0.9 * ySize, 'fwhm=%4.1f pix' % fwhmPix, ha='center', va='top')

        pltEllipseMap.set_title('Size and Ellongation of PSF sources')
        pltEllipseMap.set_xlabel('X (pix)')
        pltEllipseMap.set_ylabel('Y (pix)')
#        plt.draw()
#        plt.show()
        fname = os.path.join(qaOutputDirName, 'psfmap_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltEllipseMap


    #================================================================================================
    # making ellipticity-map plots        
    if True:
        xSize=2048.0
        ySize=4177.0
        facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
        wFig = xSize * facSize * 1.3
        hFig = ySize * facSize
        fig = plt.figure(figsize=(wFig,hFig))
        
#        pltEllMap = fig.add_subplot(1,1,1)
        pltEllMap = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height
        pltEllMap.set_xlim(0, xSize)
        pltEllMap.set_ylim(0, ySize)
#        pointSize = math.pi*(2*fwhmListPsfLikeRobust/2)**2. # 10pix=2arcsec fwhm = 10 point radius
        #pointSize = fwhmListPsfLikeRobust
#        pltEllMap.plot(xListPsfLikeRobust, yListPsfLikeRobust, markersize=pointSize, marker='o', label='PSF sample')
        #pltEllMap.scatter(xListPsfLikeRobust, yListPsfLikeRobust, s=pointSize, marker='o', label='PSF sample')
        #pltEllMap.quiver(xListPsfLikeRobust, yListPsfLikeRobust, 10*IxxListPsfLikeRobust, 10*IyyListPsfLikeRobust, headlength=0)
#        print '*** ellPaListPsfLikeRobust:', ellPaListPsfLikeRobust
#        print '*** len(ellPaListPsfLikeRobust):', len(ellPaListPsfLikeRobust)
#        for iii, xxx in enumerate(ellPaListPsfLikeRobust):
#            print iii, xxx 
#        for ell, ellPa in zip(ellListPsfLikeRobust, ellPaListPsfLikeRobust):
#            print ell, ellPa
            
#        ellPaRadianListPsfLikeRobust = numpy.array([numpy.radians(x) for x in ellPaListPsfLikeRobust if x is not None])
        ellX = numpy.array([ell*numpy.cos(numpy.radians(ellPa)) for ell, ellPa in zip(ellListPsfLikeRobust, ellPaListPsfLikeRobust) if all([ell, ellPa])])
        ellY = numpy.array([ell*numpy.sin(numpy.radians(ellPa)) for ell, ellPa in zip(ellListPsfLikeRobust, ellPaListPsfLikeRobust) if all([ell, ellPa])])
        
        #scaleFactor = 50.
        Q = pltEllMap.quiver(
            xListPsfLikeRobust, yListPsfLikeRobust,
#            ellListPsfLikeRobust*numpy.cos(ellPaRadianListPsfLikeRobust),
#            ellListPsfLikeRobust*numpy.sin(ellPaRadianListPsfLikeRobust),
            ellX, ellY,
            units = 'x', # scale is based on multiplication of 'x (pixel)'
            angles = 'uv',
            #angles = ellPaMed,
            scale = 0.0005,   # (ell/pix)
            # 1pix corresponds to this value of ellipticity --> 1000pix = ellipticity 0.02 for full CCD size.
            # and scaled by using (GridSize(in shorter axis) / CcdSize)
            scale_units = 'x',
            headwidth=0,
            headlength=0,
            headaxislength=0,
            label='PSF sample'
            )

        plt.quiverkey(Q, 0.05, 1.05, 0.05, 'e=0.05', labelpos='W')

        pltEllMap.set_title('Ellipticity of PSF sources')
        pltEllMap.set_xlabel('X (pix)')
        pltEllMap.set_ylabel('Y (pix)')
#        plt.draw()
        pltEllMap.legend()
#        plt.show()
        fname = os.path.join(qaOutputDirName, 'ellmap_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltEllMap

    #================================================================================================
    # making fwhm-in-grids plots

    xCcd = 2048.
    yCcd = 4177.
    xGridSize = yGridSize = 1024.

    outputFileName = os.path.join(qaOutputDirName, 'fwhmccd_'+frameId+'.png')
    makePlotFwhmInGrids(xCcd, yCcd, xGridSize, yGridSize, xListPsfLikeRobust, yListPsfLikeRobust, fwhmListPsfLikeRobust, outputFileName = outputFileName)

    #================================================================================================
    # making ellipse-in-grids plots

    xCcd = 2048.
    yCcd = 4177.
    xGridSize = yGridSize = 1024.

    outputFileName = os.path.join(qaOutputDirName, 'psfccd_'+frameId+'.png')
    makePlotPsfEllipseInGrids(xCcd, yCcd, xGridSize, yGridSize, xListPsfLikeRobust, yListPsfLikeRobust, AEllListPsfLikeRobust, BEllListPsfLikeRobust, ellPaListPsfLikeRobust, outputFileName = outputFileName)

    #================================================================================================
    # making ellipticity-in-grids plots

    # http://stackoverflow.com/questions/9382664/python-matplotlib-imshow-custom-tickmarks
    
    if True:
        xSize=2048.0
        ySize=4177.0
        facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
#        wFig = xSize * facSize * 1.3
        wFig = xSize * facSize * 1.2
        hFig = ySize * facSize
        fig = plt.figure(figsize=(wFig,hFig))
        
        pltEll = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height
        pltEll.set_xlim(0, xSize)
        pltEll.set_ylim(0, ySize)

        # grid Size
        if True:
            xGridSize = 1024.0
            yGridSize = 1024.0
            nx = int(numpy.floor(xSize/xGridSize))
            ny = int(numpy.floor(ySize/yGridSize))

            xGrids = numpy.array([ (i+0.5)*xGridSize for i in range(nx) ])
            yGrids = numpy.array([ (i+0.5)*yGridSize for i in range(ny) ])
            xMeshGrid, yMeshGrid = numpy.meshgrid(xGrids, yGrids)
            xGridList = numpy.reshape(xMeshGrid, (-1,))
            yGridList = numpy.reshape(yMeshGrid, (-1,))            
            
            ellXList = numpy.array([])
            ellYList = numpy.array([])
            ellMed = numpy.array([])
            ellPaMed = numpy.array([])
            # walking through all grid meshes
            for xGridCenter, yGridCenter in zip(xGridList, yGridList):
                xGridMin = xGridCenter - xGridSize/2.
                xGridMax = xGridCenter + xGridSize/2.
                yGridMin = yGridCenter - yGridSize/2.
                yGridMax = yGridCenter + yGridSize/2.

                #print '**xmin, xmax, ymin, ymax:', xGridMin, xGridMax, yGridMin, yGridMax
                ellsInGrid = numpy.array([])
                ellPasInGrid = numpy.array([])
                for xEll, yEll, ell, ellPa in zip(xListPsfLikeRobust, yListPsfLikeRobust, ellListPsfLikeRobust, ellPaListPsfLikeRobust):
                    #print '**xEll,yEll, ell, ellPa:', xEll, yEll, ell, ellPa
                    if (xGridMin <= xEll) and (xEll < xGridMax) and (yGridMin <= yEll) and (yEll < yGridMax):
                        if all([ell, ellPa]):
                            ellsInGrid = numpy.append(ellsInGrid, ell)
                            ellPasInGrid = numpy.append(ellPasInGrid, ellPa)
                # taking median to represent the value in a grid mesh
                ellPerGrid = numpy.median(ellsInGrid)
                ellPaPerGrid = numpy.median(ellPasInGrid)
#                if math.isnan(ellPerGrid) or math.isnan(ellPaPerGrid):
#                    ellPerGrid = 0.0
#                    ellPaPerGrid = 0.0                    

                ellXList = numpy.append(ellXList, ellPerGrid*numpy.cos(numpy.radians(ellPaPerGrid)))
                ellYList = numpy.append(ellYList, ellPerGrid*numpy.sin(numpy.radians(ellPaPerGrid)))

                ellMed = numpy.append(ellMed, ellPerGrid)
                ellPaMed = numpy.append(ellPaMed, ellPaPerGrid)                
            #print '** ellXList, ellYList:', ellXList, ellYList
        else:                                                                
            # testing case. only 1 grid per ccd
            xGridSize, yGridSize = xSize, ySize
            xGridList = [xSize/2.]
            yGridList = [ySize/2.]        
            ellXList = numpy.array([ ellRobust*numpy.cos(numpy.radians(ellPaRobust)) ])
            ellYList = numpy.array([ ellRobust*numpy.sin(numpy.radians(ellPaRobust)) ])
            for ellx, elly in zip(ellXList, ellYList):
                print '**** (ellX, ellY) = ( %f, %f )' % (ellx, elly)

        scaleFactor = min(xGridSize/xSize, yGridSize/ySize)
        print '**** scaleFactor', scaleFactor

        Q = pltEll.quiver(
            xGridList, yGridList,
            ellXList, ellYList,    # ellipticity 1 ~ 1000?
#            10000*ellX, 10000*ellY,    # ellipticity 1 ~ 1000?        
#            min(xGridSize, yGridSize)*0.5*ellX, min(xGridSize, yGridSize)*0.5*ellY,    # ellipticity 1 ~ 1000?
            units = 'x', # scale is based on multiplication of 'x (pixel)'
            angles = 'uv',
            #angles = ellPaMed,
            scale = (2e-5 / scaleFactor),   # (ell/pix)
            # 1pix corresponds to this value of ellipticity --> 1000pix = ellipticity 0.02 for full CCD size.
            # and scaled by using (GridSize(in shorter axis) / CcdSize)
            scale_units = 'x',
            headwidth=0,
            headlength=0,
            headaxislength=0,
            width = 50,
            #width = 0.005 * min(xGridSize, yGridSize),
#            label='PSF sample'
            #angles = 'xy',
            pivot = 'middle', #pivot = 'tail'
            )
        plt.quiverkey(Q, 0.05, 1.05, 0.05, 'e=0.05', labelpos='W')

        plt.xticks([ xc+xGridSize/2. for xc in xGridList ])
        plt.yticks([ yc+yGridSize/2. for yc in yGridList ])
        pltEll.grid()
        
        pltEll.set_title('Ellipticity of PSF sources')
        pltEll.set_xlabel('X (pix)')
        pltEll.set_ylabel('Y (pix)')
        #plt.draw()
        #pltEll.legend()
        #plt.show()
        fname = os.path.join(qaOutputDirName, 'ellccd_'+frameId+'.png')
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltEll


    #================================================================================================
    
    #print '*** fwhmRobust: %f (pix) SC: %f (arcsec) HSC: %f (arcsec)  ellRobust: %f' % (fwhmRobust, fwhmRobust*0.202, fwhmRobust*0.168, ellRobust)

    # preparing sourceSet which has been used for rough FWHM estimation
    #sourceSetPsfLike = afwDetection.SourceSet()
    #sourceSetPsfLikeRobust = afwDetection.SourceSet()

    # catalogPsfLike = afwTable.Catalog.Catalog(catalogSchema)
    ##### It would be nice if I can add measured fwhm, ellipticity, pa of each source entry here 
    catalogPsfLike = afwTable.SourceCatalog(catalogSchema)
    catalogPsfLikeRobust = afwTable.SourceCatalog(catalogSchema)
    #for i, record in enumerate(catalog):
    #    print i, record
    #sys.exit(0)

    ## note by FH: SourceCatalog object seems to not accept STL vector-like functions
    ## e.g., push_back, reserve etc, though those are listed in doxygen member functions. my misunderstanding?
    ## Instead, I use here addNew() and copy record reference to the new one.
    if False:
        catalogPsfLike.reserve(len(indicesSourcePsfLike))
        catalogPsfLikeRobust.reserve(len(indicesSourcePsfLikeRobust))        

        for i, j in enumerate(indicesSourcesPsfLike):
            catalogPsfLike[i] = catalog.table.copyRecord(catalog[j])
        for i, j in enumerate(indicesSourcesPsfLikeRobust):
            catalogPsfLikeRobust[i] = catalog.table.copyRecord(catalog[j])
    elif True:
        for i, record in enumerate(catalog):
            if i in indicesSourcesPsfLike:
                #catalogPsfLike.push_back(record)
                recordNew = catalogPsfLike.addNew()
                recordNew = record  # do we need catalog.table.copyRecord(record) here?
            if i in indicesSourcesPsfLikeRobust:
                #catalogPsfLikeRobust.push_back(record)
                recordNew = catalogPsfLikeRobust.addNew()
                recordNew = record

    # filling up the metadata keywords for QA
    metadata = exposure.getMetadata()
    metadata.set('SEEING_MODE', fwhmRobust)
    metadata.set('ELL_MED', ellRobust)
    metadata.set('ELL_PA_MED', ellPaRobust)    

    return fwhmRobust, ellRobust, ellPaRobust, catalogPsfLike, catalogPsfLikeRobust


def getVisitIdAndCcdIdFromFrameId(frameId, config):
    span = re.search('[0-9]*[0-9]', frameId).span() # extracting a numeric part
    if False:
        camera = 'suprimecam'
    else:
        #camera = config['camera']
        #camera = config.camera
        camera = config.qa.camera # I need know here what the instrument is or need any other way to get visitId and ccdId

    print '*** measSeeing: camera:', camera
    if camera.lower() == "hsc" or camera.lower() == "hscsim":
        visitId = int(frameId[span[0]:span[1]-3])
        ccdId = int(frameId[span[1]-3:])
        print '*** measSeeingQa:getVisitIdAndCcdIdFromFrameId - Hsc camera: %s visitId: %d  ccdId:  %d' % (camera, visitId, ccdId)        
    elif camera.lower() in ("suprimecam", "suprime-cam", "sc", "suprimecam-mit", "sc-mit", "scmit", "suprimecam-old", "sc-old", "scold"):
        visitId = int(frameId[span[0]:span[1]-1])
        ccdId = int(frameId[-1])
        print '*** measSeeingQa:getVisitIdAndCcdIdFromFrameId - Sc camera: %s visitId: %d  ccdId:  %d' % (camera, visitId, ccdId)
    else:
        print 'Instrument specified is invalid.'

    return visitId, ccdId


def makePlotFwhmInGrids(xCcd, yCcd, xGridSize, yGridSize, xListPsfLikeRobust, yListPsfLikeRobust, fwhmListPsfLikeRobust, outputFileName = None):
    xSize = xCcd
    ySize = yCcd
    facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
    wFig = xSize * facSize * 1.3
    hFig = ySize * facSize
    fig = plt.figure(figsize=(wFig,hFig))
    pltFwhm = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height

    pltFwhm.set_xlim(0, xSize)
    pltFwhm.set_ylim(0, ySize)

    # making grids
    nx = int(numpy.floor(xSize/xGridSize))
    ny = int(numpy.floor(ySize/yGridSize))

    xGrids = numpy.array([ (i+0.5)*xGridSize for i in range(nx) ])
    yGrids = numpy.array([ (i+0.5)*yGridSize for i in range(ny) ])
    xMeshGrid, yMeshGrid = numpy.meshgrid(xGrids, yGrids)
    xGridList = numpy.reshape(xMeshGrid, (-1,))
    yGridList = numpy.reshape(yMeshGrid, (-1,))    

    # walking through all grid meshes
    fwhmList = numpy.array([])
    for xGridCenter, yGridCenter in zip(xGridList, yGridList):
        xGridMin = xGridCenter - xGridSize/2.
        xGridMax = xGridCenter + xGridSize/2.
        yGridMin = yGridCenter - yGridSize/2.
        yGridMax = yGridCenter + yGridSize/2.

        #print '**xmin, xmax, ymin, ymax:', xGridMin, xGridMax, yGridMin, yGridMax
        fwhmsInGrid = numpy.array([])
        for xFwhm, yFwhm, fwhm in zip(xListPsfLikeRobust, yListPsfLikeRobust, fwhmListPsfLikeRobust):
            #print '**xFwhm, yFwhm, fwhm:', xFwhm, yFwhm, fwhm
            if (xGridMin <= xFwhm) and (xFwhm < xGridMax) and (yGridMin <= yFwhm) and (yFwhm < yGridMax):
                if fwhm is not None:
                    fwhmsInGrid = numpy.append(fwhmsInGrid, fwhm)
        # taking median to represent the value in a grid mesh
        fwhmList = numpy.append(fwhmList, numpy.median(fwhmsInGrid))

    print '** fwhmList:', fwhmList

    pointRadius = 100*fwhmList/2. # 10pix=2arcsec(fwhm)=500 point(radius) (to be ~0.6*min(xGridSize, yGridSize)?)
    scaleFactor = min(xGridSize/xSize, yGridSize/ySize)
    pointRadius *= scaleFactor    
    pointArea = math.pi*(pointRadius)**2.

    #pltFwhm.plot(xListPsfLikeRobust, yListPsfLikeRobust, markersize=pointSize, marker='o', label='PSF sample')
    pltFwhm.scatter(xGridList, yGridList, s=pointArea, marker='o', color=None, facecolor=(1,1,1,0), linewidth=5.0, label='PSF sample')

    # reference sample symbol
    fwhmPix = 5.0 # 5pix in fwhm
    pointRadius = 100*numpy.array([fwhmPix])/2. 
    scaleFactor = min(xGridSize/xSize, yGridSize/ySize)
    pointRadius *= scaleFactor    
    pointArea = math.pi*(pointRadius)**2.
    pltFwhm.scatter([0.1 * xSize], [0.9 * ySize], s=pointArea, marker='o', color='magenta', facecolor=(1,1,1,0), linewidth=8.0, label='PSF sample')
    plt.text(0.1 * xSize, 0.9 * ySize, 'fwhm=%4.1f pix' % fwhmPix, ha='center', va='top')
    
    pltFwhm.set_title('FWHM of PSF sources')
    pltFwhm.set_xlabel('X (pix)')
    pltFwhm.set_ylabel('Y (pix)')
    #plt.draw()
    #pltFwhm.legend()

    plt.xticks([ xc+xGridSize/2. for xc in xGridList ])
    plt.yticks([ yc+yGridSize/2. for yc in yGridList ])
    pltFwhm.grid()
    
    #plt.show()
    if not outputFileName:
        outputFileName = ''
    plt.savefig(outputFileName, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

    del fig
    del pltFwhm

def makePlotPsfEllipseInGrids(xCcd, yCcd, xGridSize, yGridSize, xListPsfLikeRobust, yListPsfLikeRobust, AEllListPsfLikeRobust, BEllListPsfLikeRobust, ellPaListPsfLikeRobust, outputFileName = None):
    xSize = xCcd
    ySize = yCcd
    facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
    wFig = xSize * facSize * 1.2
    hFig = ySize * facSize
    fig = plt.figure(figsize=(wFig,hFig))
    pltEllipse = plt.axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height

    pltEllipse.set_xlim(0, xSize)
    pltEllipse.set_ylim(0, ySize)

    # making grids
    nx = int(numpy.floor(xSize/xGridSize))
    ny = int(numpy.floor(ySize/yGridSize))

    xGrids = numpy.array([ (i+0.5)*xGridSize for i in range(nx) ])
    yGrids = numpy.array([ (i+0.5)*yGridSize for i in range(ny) ])
    xMeshGrid, yMeshGrid = numpy.meshgrid(xGrids, yGrids)
    xGridList = numpy.reshape(xMeshGrid, (-1,))
    yGridList = numpy.reshape(yMeshGrid, (-1,))    

    # walking through all grid meshes
    AEllList = numpy.array([])
    BEllList = numpy.array([])
    ellPaList = numpy.array([])        
    for xGridCenter, yGridCenter in zip(xGridList, yGridList):
        xGridMin = xGridCenter - xGridSize/2.
        xGridMax = xGridCenter + xGridSize/2.
        yGridMin = yGridCenter - yGridSize/2.
        yGridMax = yGridCenter + yGridSize/2.

        #print '**xmin, xmax, ymin, ymax:', xGridMin, xGridMax, yGridMin, yGridMax
        AEllsInGrid = numpy.array([])
        BEllsInGrid = numpy.array([])        
        ellPasInGrid = numpy.array([])
        for xEllipse, yEllipse, aa, bb, ellPa in zip(xListPsfLikeRobust, yListPsfLikeRobust, AEllListPsfLikeRobust, BEllListPsfLikeRobust, ellPaListPsfLikeRobust):
            #print '**xPsf, yPsf, fwhm:', xPsf, yPsf, fwhm
            if (xGridMin <= xEllipse) and (xEllipse < xGridMax) and (yGridMin <= yEllipse) and (yEllipse < yGridMax):
                if all([aa, bb, ellPa]):
                    AEllsInGrid = numpy.append(AEllsInGrid, aa)
                    BEllsInGrid = numpy.append(BEllsInGrid, bb)
                    ellPasInGrid = numpy.append(ellPasInGrid, ellPa)
        # taking median to represent the value in a grid mesh
        AEllList = numpy.append(AEllList, numpy.median(AEllsInGrid))
        BEllList = numpy.append(BEllList, numpy.median(BEllsInGrid))        
        ellPaList = numpy.append(ellPaList, numpy.median(ellPasInGrid))        
    print '** AEllList:', AEllList
    print '** BEllList:', BEllList    

    #pltPsf.plot(xListPsfLikeRobust, yListPsfLikeRobust, markersize=pointSize, marker='o', label='PSF sample')
    scaleFactor = (1/10.)*0.8*min(xGridSize, yGridSize) # 10pix=2arcsec(fwhm)==>A=0.6*gridSize~1200pix(for whole_CCD)
    for xEllipse, yEllipse, aa, bb, ellPa in zip(xGridList, yGridList, AEllList, BEllList, ellPaList):
        if all([aa, bb, ellPa]):
            ell = patches.Ellipse((xEllipse, yEllipse), 2.*aa*scaleFactor, 2.*bb*scaleFactor, angle=ellPa, linewidth=2., fill=False, zorder=2)
        pltEllipse.add_patch(ell)
    pltEllipse.set_title('Size and Ellongation of PSF sources')
    pltEllipse.set_xlabel('X (pix)')
    pltEllipse.set_ylabel('Y (pix)')
    #plt.draw()
    #pltEllipse.legend()

    # reference sample point
    fwhmPix = 5 # pixel in A, B (in half width)
    aa = bb = fwhmPix/2.
    ell = patches.Ellipse((0.1*xSize, 0.9*ySize), 2.*aa*scaleFactor, 2.*bb*scaleFactor, angle=0., linewidth=4., color='magenta', fill=False, zorder=2)
    pltEllipse.add_patch(ell)
    plt.text(0.1 * xSize, 0.9 * ySize, 'fwhm=%4.1f pix' % fwhmPix, ha='center', va='top')

    plt.xticks([ xc+xGridSize/2. for xc in xGridList ])
    plt.yticks([ yc+yGridSize/2. for yc in yGridList ])
    pltEllipse.grid()
    
    #plt.show()
    if not outputFileName:
        outputFileName = ''
    plt.savefig(outputFileName, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

    del fig
    del pltEllipse

