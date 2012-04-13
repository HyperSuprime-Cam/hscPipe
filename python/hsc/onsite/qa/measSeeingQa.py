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
plt.switch_backend('Agg') # experimental option to work around unknown conflict in the Agg backend 
### Be sure that other python modules do not call matplotlib.pyplot with TkAgg or non-Agg backend before this module.


class QaSeeingConfig(pexConfig.Config):
    fwhmIni = pexConfig.Field(
        dtype = float,
        doc = 'Inintial guess of seeing (pix)',
        default = 3.465,
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
    
def measureSeeingQa(exposure, catalog, config, debugFlag=False, plotFlag=True, plotbasedir=None, butler=None, log=None):
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


    # -- Filtering sources with rough min-max fwhm
    magListAll = []
    fwhmListAll = []
    ellListAll = []
    ellPaListAll = []        
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

        Ixx = source.getIxx() # luminosity-weighted 2nd moment of pixels
        Iyy = source.getIyy() 
        Ixy = source.getIxy()
        sigma = math.sqrt(0.5*(Ixx+Iyy))
        fwhm = 2.*math.sqrt(2.*math.log(2.)) * sigma # (pix) assuming Gaussian

        saturationFlag = source.get(keySaturationFlag)
        saturationCenterFlag = source.get(keySaturationCenterFlag)        
#        objFlag = source.getFlagForDetection()
        isSaturated = (saturationFlag | saturationCenterFlag)        
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

        print '*** %d : Ixx: %f Iyy: %f Ixy: %f fwhm: %f flux: %f isSatur: %x' % (iseq, Ixx, Iyy, Ixy, fwhm, fluxAper, isSaturated)
        mag = -2.5*numpy.log10(fluxForSeeing)
        magListAll.append(mag)
        fwhmListAll.append(fwhm)
        
        if True: # definition by SExtractor
            Val1 = 0.5*(Ixx+Iyy)
            Ixx_Iyy = Ixx-Iyy
            Val2 = 0.25*Ixx_Iyy*Ixx_Iyy + Ixy*Ixy
            if Val2 >= 0 and (Val1-math.sqrt(Val2)) > 0 and Ixx_Iyy != 0:
                aa = math.sqrt( Val1 + math.sqrt(Val2) )
                bb = math.sqrt( Val1 - math.sqrt(Val2) )
                ell =  1. - bb/aa
                ellPa = 0.5 * math.degrees(math.atan(2*Ixy / math.fabs(Ixx_Iyy)))
            else:
                ell = None
                ellPa = None
        else: # definition by Kaiser
            # e=sqrt(e1^2+e2^2) where e1=(Ixx-Iyy)/(Ixx+Iyy), e2=2Ixy/(Ixx+Iy)
            # SExtractor's B/A=sqrt((1-e)/(1+e)), ell=1-B/A
            e1 = (Ixx-Iyy)/(Ixx+Iyy)
            if e1 > 0: 
                e2 = 2.0*Ixy/(Ixx+Iyy)
                ell = math.sqrt(e1*e1 + e2*e2)
                ellPa = 0.5 * math.degrees(math.atan(2*Ixy / math.fabs(Ixx-Iyy)))            
            else:
                ell = 0
                ellPa = 0

        ellListAll.append( ell )
        ellPaListAll.append( ellPa )
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
        if cumFraction >= 0.15:
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

    # making debug plots
    if plotFlag is True: 
        log.log(log.INFO, "QaSeeing: Entered plotting")
        # getting visitId and ccdId. here, we reconstruct those two values from FrameId.
        # maybe we should consider a better way.
        frameId = exposure.getMetadata().get('FRAMEID')
        #print 'XXXXXXX QaSeeing instrument = ', instrument
        visitId, ccdId = getVisitIdAndCcdIdFromFrameId(frameId, config)

        if False:
            det = exposure.getDetector()
            ccd = cameraGeom.cast_Ccd(det)
            ccdId = int(ccd.getId().getSerial())

        log.log(log.INFO, "QaSeeing: After getting ccdId")

        dataId = { 'visit': visitId, 'ccd': ccdId }
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

        fname = qaOutputDirName + '/' + 'seeing1_'+frameId+'.png'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

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
    
    #print '*** indicesSourcesFwhmRange:', indicesSourcesFwhmRange
    #print '*** indicesSoucresPsfLike:',  indicesSourcesPsfLike
    #print '*** magListPsfLike:',  magListPsfLike
    #print '*** fwhmListPsfLike:',  fwhmListPsfLike

    log.log(log.INFO, "QaSeeing: nSampleRoughFwhm: %d" % nSampleRoughFwhm)

    fwhmListForRoughFwhm = numpy.sort(fwhmListPsfLike)[:nSampleRoughFwhm] # fwhmListPsfLike is not changed with numpy.sort()
    fwhmRough = numpy.median(fwhmListForRoughFwhm)
    print '*** fwhmListPsfLike:', fwhmListPsfLike
    print '*** fwhmListForRoughFwhm:', fwhmListForRoughFwhm
    print '*** fwhmRough:', fwhmRough
    log.log(log.INFO, "QaSeeing: fwhmRough: %f" % fwhmRough)
    
    # making debug plots
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
        fname = qaOutputDirName + '/' + 'seeing2_'+frameId+'.png'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

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
        fwhmRobust = None
        ellRobust = None
        ellPaRobust = None
    else:
        magListPsfLikeRobust = magListAll[indicesSourcesPsfLikeRobust]
        fwhmListPsfLikeRobust = fwhmListAll[indicesSourcesPsfLikeRobust]
        ellListPsfLikeRobust = ellListAll[indicesSourcesPsfLikeRobust]
        ellPaListPsfLikeRobust = ellPaListAll[indicesSourcesPsfLikeRobust]        
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

    # making debug plots
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
        fname = qaOutputDirName + '/' + 'seeing3_'+frameId+'.png'
        plt.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

    #print '*** fwhmRobust: %f (pix) SC: %f (arcsec) HSC: %f (arcsec)  ellRobust: %f' % (fwhmRobust, fwhmRobust*0.202, fwhmRobust*0.168, ellRobust)

    # preparing sourceSet which has been used for rough FWHM estimation
    #sourceSetPsfLike = afwDetection.SourceSet()
    #sourceSetPsfLikeRobust = afwDetection.SourceSet()

    # catalogPsfLike = afwTable.Catalog.Catalog(catalogSchema)
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
    if camera.lower() in ("hsc", "hscsim"):
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

