#!/usr/bin/env python

import os
import lsst.pex.config as pexConfig
import lsst.afw.detection as afwDet
import lsst.afw.table as afwTable
import lsst.daf.base as dafBase
import lsst.meas.algorithms as measAlg
import lsst.pipe.base as pipeBase
import lsst.pipe.tasks.processCcd as ptProcessCcd
import hsc.pipe.tasks.astrometry as hscAstrom
import hsc.pipe.tasks.suprimecam as hscSuprimeCam
import hsc.pipe.tasks.calibrate as hscCalibrate
import hsc.pipe.tasks.isr as hscIsr
import hsc.pipe.tasks.hscDc2 as hscDc2

##== FH added for QA output
import hsc.pipe.tasks.qaSuprimeCamIsrTask as qaSuprimeCamIsrTask
import hsc.onsite.qa.measSeeingQa as qaSeeing
import hsc.pipe.tasks.qaCalibrate as qaHscCalibrate
import lsst.afw.geom as afwGeom
import numpy

#== FH changed for QA output

#class QaWriteFits(pexConfig):
class QaConfig(pexConfig.Config):
    seeing = pexConfig.ConfigField(dtype=qaSeeing.QaSeeingConfig, doc="Qa.measSeeing")
    camera = pexConfig.Field(dtype=str, doc="camera type [hsc] or [suprimecam]", default='hsc') # would be better to have another way to get camera

class SubaruProcessCcdConfig(ptProcessCcd.ProcessCcdConfig):
    calibrate = pexConfig.ConfigField(dtype=qaHscCalibrate.HscCalibrateConfig, doc="Calibration")
##    calibrate = pexConfig.ConfigField(dtype=hscCalibrate.HscCalibrateConfig, doc="Calibration")    
    isr = pexConfig.ConfigField(dtype=qaSuprimeCamIsrTask.QaIsrTaskConfig, doc="Config for Isr with Qa tasks")
    qa = pexConfig.ConfigField(dtype=QaConfig, doc="Config for Qa (outside isr) tasks")

class SubaruProcessCcdTask(ptProcessCcd.ProcessCcdTask):
    """Subaru version of ProcessCcdTask, with method to write outputs
    after producing a new multi-frame WCS.
    """
    ConfigClass = SubaruProcessCcdConfig

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdTask, self).__init__(*args, **kwargs)
        self.makeSubtask("isr", hscIsr.HscIsrTask)
        self.makeSubtask("calibrate", hscCalibrate.HscCalibrateTask)
        self.makeSubtask("photometry", PhotometryTask) # this class seems to be not exist in this script, yet

    # The 'run' method is copied wholesale from lsst.pipe.tasks.processCcd.ProcessCcdTask.run, with minor
    # modifications to change when the CCD assembly is performed.
    @pipeBase.timeMethod
    def run(self, sensorRef):
        self.log.log(self.log.INFO, "Processing %s" % (sensorRef.dataId))
        if self.config.doIsr:
            butler = sensorRef.butlerSubset.butler
            calibSet = self.isr.makeCalibDict(butler, sensorRef.dataId)
            exposure = sensorRef.get("raw")
            isrRes = self.isr.run(exposure, calibSet)
            exposure = isrRes.postIsrExposure
            self.display("isr", exposure=exposure, pause=True)
            if self.config.doWriteIsr:
                sensorRef.put(exposure, 'postISRCCD')
        else:
            exposure = None

        if self.config.doCalibrate:
            if exposure is None:
                exposure = sensorRef.get('postISRCCD')
            calib = self.calibrate.run(exposure)
            exposure = calib.exposure
            if self.config.doWriteCalibrate:
                sensorRef.put(exposure, 'calexp')
                sensorRef.put(calib.sources, 'icSrc')
                if calib.psf is not None:
                    sensorRef.put(calib.psf, 'psf')
                if calib.apCorr is not None:
                    sensorRef.put(calib.apCorr, 'apCorr')
                if calib.matches is not None:
                    normalizedMatches = afwTable.packMatches(calib.matches)
                    normalizedMatches.table.setMetadata(calib.matchMeta)
                    sensorRef.put(normalizedMatches, 'icMatch')
        else:
            calib = None

        if self.config.doDetection:
            if exposure is None:
                exposure = sensorRef.get('calexp')
            if calib is None:
                psf = sensorRef.get('psf')
                exposure.setPsf(sensorRef.get('psf'))
            table = afwTable.SourceTable.make(self.schema)
            table.setMetadata(self.algMetadata)
            detRet = self.detection.makeSourceCatalog(table, exposure)
            sources = detRet.sources
        else:
            sources = None

        if self.config.doMeasurement:
            assert(sources)
            assert(exposure)
            if calib is None:
                apCorr = sensorRef.get("apCorr")
            else:
                apCorr = calib.apCorr
            self.measurement.run(exposure, sources, apCorr)
 
        if self.config.doWriteSources:
            sensorRef.put(sources, 'src')

        if self.config.doWriteCalibrate:
            sensorRef.put(exposure, 'calexp')
            
        return pipeBase.Struct(
            exposure = exposure,
            calib = calib,
            sources = sources,
        )

    def write(self, butler, dataId, struct, wcs=None):

        if wcs is None:
            wcs = struct.exposure.getWcs()
            self.log.log(self.log.WARN, "WARNING: No new WCS provided")

        # Apply WCS to sources
        # No longer handling matchSources explicitly - these should all be in calib.sources,
        # or there's a bug in the calibrate task.
        struct.exposure.setWcs(wcs)
        for sources in (struct.sources, struct.calib.sources):
            if sources is None:
                continue
            for s in sources:
                s.updateCoord(wcs)
                
        normalizedMatches = afwTable.packMatches(struct.calib.matches)
        normalizedMatches.table.setMetadata(struct.calib.matchMeta)

        butler.put(struct.exposure, 'calexp', dataId)
        butler.put(struct.sources, 'src', dataId)
        butler.put(normalizedMatches, 'icMatch', dataId)
        butler.put(struct.calib.psf, 'psf', dataId)
        butler.put(struct.calib.apCorr, 'apCorr', dataId)
        butler.put(struct.calib.sources, 'icSrc', dataId)


class SuprimeCamProcessCcdTask(SubaruProcessCcdTask):
    
    def __init__(self, **kwargs):
        ## I would not like to call IsrTask.__init__() where methodList is created
        pipeBase.Task.__init__(self, **kwargs) 
        ## FH changed for QA outputs
        ##        self.makeSubtask("isr", hscSuprimeCam.SuprimeCamIsrTask)
        ##        self.makeSubtask("calibrate", hscCalibrate.HscCalibrateTask)
        self.makeSubtask("isr", qaSuprimeCamIsrTask.QaSuprimeCamIsrTask)#, config=SubaruProcessCcdConfig())
        self.makeSubtask("calibrate", qaHscCalibrate.HscCalibrateTask) #, config=self.config)        
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.algMetadata = dafBase.PropertyList()
        if self.config.doDetection:
            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
        if self.config.doMeasurement:
            self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
                             schema=self.schema, algMetadata=self.algMetadata)

    ## FH added the below run() function for QA output
    ## This is alomst same as in the SubaruProcessCcdTask.run() but 
    ## overriding it with minor modifications to change arguments of
    ## self.isr.run() to accept butler and dataId.
    @pipeBase.timeMethod
    def run(self, sensorRef):
        self.log.log(self.log.INFO, "Processing %s" % (sensorRef.dataId))
        if self.config.doIsr:
            butler = sensorRef.butlerSubset.butler
            calibSet = self.isr.makeCalibDict(butler, sensorRef.dataId)
            exposure = sensorRef.get("raw")
            dataId = sensorRef.dataId
            ##isrRes = self.isr.run(rawExposure, calibSet)
            ##isrRes = self.isr.run(rawExposure, calibSet, butler, dataId)            
            isrRes = self.isr.run(exposure, calibSet, butler, dataId)            
            ##exposure = self.isr.doCcdAssembly([isrRes.postIsrExposure])
            exposure = isrRes.postIsrExposure
            self.display("isr", exposure=exposure, pause=True)
            if self.config.doWriteIsr:
                sensorRef.put(exposure, 'postISRCCD')
        else:
            exposure = None

        if self.config.doCalibrate:
            if exposure is None:
                exposure = sensorRef.get('postISRCCD')
            calib = self.calibrate.run(exposure)
            exposure = calib.exposure
            if self.config.doWriteCalibrate:
                sensorRef.put(exposure, 'calexp')
                sensorRef.put(calib.sources, 'icSrc')
                if calib.psf is not None:
                    sensorRef.put(calib.psf, 'psf')
                if calib.apCorr is not None:
                    sensorRef.put(calib.apCorr, 'apCorr')
                if calib.matches is not None:
                    normalizedMatches = afwTable.packMatches(calib.matches)
                    print '*** len(normalizedMatches):',len(normalizedMatches)
                    normalizedMatches.table.setMetadata(calib.matchMeta)
                    sensorRef.put(normalizedMatches, 'icMatch')

                    # === For QA db
                    sdssFluxes = retrieveReferenceMagnitudes(calib.matches)
                    #print '*** returned sdssFluxes:', sdssFluxes
                    #writeMatchesToBintableFits(calib.matches, butler=sensorRef)
                    writeEnhancedMatchesToBintableFits(calib.matches, calib.matchMeta, refMags=sdssFluxes, butler=sensorRef)                    
        else:
            calib = None

        if True:
            ##== FH added this part for QA output to use the reduced exposure and (bright) sources in calibration 
            butler = sensorRef.butlerSubset.butler            
            fwhmRobust, ellRobust, ellPaRobust, catalogPsfLike, catalogPsfLikeRobust = qaSeeing.measureSeeingQa(exposure, calib.sources, self.config, debugFlag=False, plotFlag=True, plotbasedir=None, butler=butler, log=self.log)

            self.log.log(self.log.INFO, "QA seeing: fwhm: %f (pix)" % fwhmRobust)
            self.log.log(self.log.INFO, "QA seeing: ell (based on 2nd moments): %f" % ellRobust)
            self.log.log(self.log.INFO, "QA seeing: ellPa (in CCDCoords based on 2nd moments): %f (deg)" % ellPaRobust)
            self.log.log(self.log.INFO, "QA seeing: final Nsources for seeing: %d" % len(catalogPsfLikeRobust))        

        ## Final source detection
        if self.config.doDetection:
            if exposure is None:
                exposure = sensorRef.get('calexp')
            if calib is None:
                psf = sensorRef.get('psf')
                exposure.setPsf(sensorRef.get('psf'))
            table = afwTable.SourceTable.make(self.schema)
            table.setMetadata(self.algMetadata)
            detRet = self.detection.makeSourceCatalog(table, exposure)
            sources = detRet.sources
        else:
            sources = None

        if self.config.doMeasurement:
            assert(sources)
            assert(exposure)
            if calib is None:
                apCorr = sensorRef.get("apCorr")
            else:
                apCorr = calib.apCorr
            self.measurement.run(exposure, sources, apCorr)

        if False:
            ##== QA output is now to use (bright) sources in calibration rather than
            ## final sources by deep detection. So, commented measSeeing here in the meantime.
            butler = sensorRef.butlerSubset.butler            
            fwhmRobust, ellRobust, ellPaRobust, catalogPsfLike, catalogPsfLikeRobust = qaSeeing.measureSeeingQa(exposure, sources, self.config, debugFlag=False, plotFlag=True, plotbasedir=None, butler=butler, log=self.log)

            self.log.log(self.log.INFO, "QA seeing: fwhm: %f (pix)" % fwhmRobust)
            self.log.log(self.log.INFO, "QA seeing: ell (based on 2nd moments): %f" % ellRobust)
            self.log.log(self.log.INFO, "QA seeing: ellPa (in CCDCoords based on 2nd moments): %f (deg)" % ellPaRobust)
            self.log.log(self.log.INFO, "QA seeing: final Nsources for seeing: %d" % len(catalogPsfLikeRobust))        

        ## == Putting in necessary information for QA data management

        metadata = exposure.getMetadata()
        
        # = flags
        # this part should be done by calculating merit functions somewhere else in a polite manner.
        metadata.set('FLAG_AUTO', 0)
        metadata.set('FLAG_USR', 0)
        metadata.set('FLAG_TAG', 1)

        # = info for data management
        # = rerun; assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'

        rerunName = getRerunName(sensorRef)
        metadata.set('RERUN', rerunName)

        if self.config.doWriteSources and sources is not None:
        ##if self.config.doWriteSources:
            sensorRef.put(sources, 'src')

        if self.config.doWriteCalibrate:
            sensorRef.put(exposure, 'calexp')

        return pipeBase.Struct(
            exposure = exposure,
            calib = calib,
            sources = sources,
        )


def getRerunName(sensorRef):
    corrPath = sensorRef.get('calexp_filename')[0]
    rerunName = corrPath[corrPath.find('rerun'):].split('/')[1]        
    #print '*** corrPath:', corrPath
    return rerunName

#
class HscProcessCcdTask(SuprimeCamProcessCcdTask):
####class HscProcessCcdTask(SubaruProcessCcdTask):
    pass

#    def __init__(self, **kwargs):
#        pipeBase.Task.__init__(self, **kwargs)
#        self.makeSubtask("isr", hscIsr.HscIsrTask)
#        self.makeSubtask("calibrate", hscCalibrate.HscCalibrateTask)
#        self.schema = afwTable.SourceTable.makeMinimalSchema()
#        self.algMetadata = dafBase.PropertyList()
#        if self.config.doDetection:
#            self.makeSubtask("detection", measAlg.SourceDetectionTask, schema=self.schema)
#        if self.config.doMeasurement:
#            self.makeSubtask("measurement", measAlg.SourceMeasurementTask,
#                             schema=self.schema, algMetadata=self.algMetadata)


def getSolverFromAstrometryNet():
    #solver = self._getSolver()
    import lsst.meas.astrom.astrometry_net as an
    solver = an.solver_new()
    # HACK, set huge default pixel scale range.
    lo,hi = 0.01, 3600.
    solver.setPixelScaleRange(lo, hi)
    return solver

def retrieveReferenceMagnitudes(matchlist):
    '''
    Searches for reference-catalog sources (in the
    astrometry_net_data files) in the requested RA,Dec region
    (afwGeom::Angle objects), with the requested radius (also an
    Angle).  The flux values will be set based on the requested
    filter (None => default filter).
    Returns: an lsst.afw.table.SimpleCatalog of reference objects
    '''
    # instantiation of andConfig
    from lsst.meas.astrom.config import MeasAstromConfig, AstrometryNetDataConfig
    import hsc.meas.astrom.astrom as hscAstrom

    config = hscAstrom.TaburAstrometryConfig()
    astrometer = hscAstrom.TaburAstrometry(config, andConfig=None)
    
    solver = getSolverFromAstrometryNet()

    filterName = 'g' 
    sgCol = astrometer.andConfig.starGalaxyColumn
    varCol = astrometer.andConfig.variableColumn
    idcolumn = astrometer.andConfig.idColumn

    sdssFilters = ['u', 'g', 'r', 'i', 'z']
    sdssMagCols = [ astrometer.getCatalogFilterName(sdssFilter) for sdssFilter in sdssFilters ]
    sdssMagErrCols = [ astrometer.andConfig.magErrorColumnMap.get(sdssFilter, None) for sdssFilter in sdssFilters ]

    if False:
        print '*** filterName:', filterName
        print '*** magCols:', sdssMagCols
        print '*** magerrCols:',  sdssMagErrCols
        print '*** idCol:', idcolumn

    # == retrieving magnitudes of each match...
    refSchema = matchlist[0].first.schema
    #refSchema: Schema(
    #    (Field['I8'](name="id", doc="unique ID"), Key<I8>(offset=0, nElements=1)),
    #        (Field['Coord'](name="coord", doc="position in ra/dec", units="IRCS; radians"), Key<Coord>(offset=8, nElements=2)),
    #        (Field['F8'](name="flux", doc="flux"), Key<F8>(offset=24, nElements=1)),
    #        (Field['F8'](name="flux.err", doc="flux uncertainty"), Key<F8>(offset=32, nElements=1)),
    #        (Field['Flag'](name="stargal", doc="set if the reference object is a star"), Key['Flag'](offset=40, bit=0)),
    #        (Field['Flag'](name="photometric", doc="set if the reference object can be used in photometric calibration"), Key['Flag'](offset=40, bit=1)),
     #   )
    fluxKey = refSchema.find('flux').key
    fluxErrKey = refSchema.find('flux.err').key    

    sdssFluxesForMatches  = []
    for match in matchlist:

        refId = match.first.getId()
        ra = match.first.getRa().asDegrees() # float in degrees
        dec = match.first.getDec().asDegrees() # float in degrees
        radius = (1.0 / 3600.0) # float in degrees

        sdssFluxes = []
        sdssFluxErrors = []
        # do we have some API function like solver.getCatalogue & getTanAlong in the previous version?
        # Instead, I am repeating getCatalog with the cone search for one band information.
        # This part would need improvement
        for sdssMagCol, sdssMagErrCol in zip(sdssMagCols, sdssMagErrCols):
            cat = solver.getCatalog(astrometer.inds, ra, dec, radius, idcolumn, sdssMagCol, sdssMagErrCol, sgCol, varCol, True)
            for s in cat:
                catId = s.getId()
                if refId != catId:
                    continue
                sdssFluxes.append(s.get(fluxKey))
                sdssFluxErrors.append(s.get(fluxErrKey))                

        if False: # debug code
            print refId, ra, dec, match.first.get(fluxKey), match.second.getId(), match.second.getX(), match.second.getY(), -2.5*numpy.log10(match.second.getPsfFlux()), '<==>', catId, sdssFluxes, sdssFluxErrors

        sdssFluxesForMatches.append( [sdssFilters, sdssFluxes, sdssFluxErrors] ) 

    del solver
    return sdssFluxesForMatches
    
def namedCopy(dstRecord, dstName, srcRecord, srcName):
    dstKey = dstRecord.schema.find(dstName).key
    srcKey = srcRecord.schema.find(srcName).key
    dstRecord.set(dstKey, srcRecord.get(srcKey))

def writeEnhancedMatchesToBintableFits(matchlist, matchMeta, refMags=None, butler=None, fileName=None):

    # - creating a new catalog schema for output involving
    #  both ref(plus reference magnitudes in additional filters) and src
    # - writing the new catalog into a bintable FITS file

    if refMags is None:
        print '** additional reference magnitudes are not provided.'

    refSchema = matchlist[0].first.getSchema()
    srcSchema = matchlist[0].second.getSchema()
    if False:
        print 'refSchema:', refSchema.getNames()
        print 'srcSchema:', srcSchema.getNames()
    
    mergedSchema = afwTable.Schema()

    for keyName in refSchema.getNames():
        field = refSchema.find(keyName).field
        typeStr = field.getTypeString()
        fieldDoc = field.getDoc()
        fieldUnits = field.getUnits()
        mergedSchema.addField('ref.'+keyName, type=typeStr, doc=fieldDoc, units=fieldUnits)

    filters, fluxes, fluxerrs = refMags[0]
    for filter in filters:
        #print '** filter', filter
        mergedSchema.addField('ref.flux.'+filter, type='F8', doc='reference flux in %s band' % filter)
        mergedSchema.addField('ref.flux.err.'+filter, type='F8', doc='reference flux error in %s band' % filter)

    for keyName in srcSchema.getNames():
        field = srcSchema.find(keyName).field
        typeStr = field.getTypeString()
        fieldDoc = field.getDoc()
        fieldUnits = field.getUnits()
        mergedSchema.addField('src.'+keyName, type=typeStr, doc=fieldDoc, units=fieldUnits)
                
    mergedCatalog = afwTable.BaseCatalog(mergedSchema)

    refKeys = []
    for keyName in refSchema.getNames():
        refKeys.append((refSchema.find(keyName).key, mergedSchema.find('ref.' + keyName).key))
    srcKeys = []
    for keyName in srcSchema.getNames():
        srcKeys.append((srcSchema.find(keyName).key, mergedSchema.find('src.' + keyName).key))

    for match, enhancedRefMags in zip(matchlist, refMags):
        record = mergedCatalog.addNew()
        for key in refKeys:
            keyIn = key[0]
            keyOut = key[1]
            record.set(keyOut, match.first.get(keyIn))
            
        for key in srcKeys:
            keyIn = key[0]
            keyOut = key[1]
            record.set(keyOut, match.second.get(keyIn))

        # this is for reference magnitudes in additional bands
        #print '*** enhancedRefMags:', enhancedRefMags
        filters, fluxes, fluxerrs = enhancedRefMags
        #print '*** filters, fluxes, fluxerrs:', filters, fluxes, fluxerrs
        for filter, flux, fluxerr in zip(filters, fluxes, fluxerrs):
            record.set(mergedSchema.find('ref.flux.'+filter).key, flux)
            record.set(mergedSchema.find('ref.flux.err.'+filter).key, fluxerr)            

    # obtaining reference catalog name
    catalogName = os.path.basename(os.getenv("ASTROMETRY_NET_DATA_DIR").rstrip('/'))
    matchMeta.add('REFCAT', catalogName)
    mergedCatalog.getTable().setMetadata(matchMeta)

    if butler is not None:
        butler.put(mergedCatalog, 'matchedList')
        return 
    elif fileName is not None:
        mergedCatalog.writeFits(fileName)
        return
    else:
        print '** No file has been written.'
        return

def writeMatchesToBintableFits(matchlist, butler=None, fileName=None):

    refSchema = matchlist[0].first.getSchema()
    srcSchema = matchlist[0].second.getSchema()
    if False:
        print 'refSchema:', refSchema.getNames()
        print 'srcSchema:', srcSchema.getNames()

    # creating a new catalog schema for output involving both ref and src
    mergedSchema = afwTable.Schema()

    for keyName in refSchema.getNames():
        field = refSchema.find(keyName).field
        typeStr = field.getTypeString()
        mergedSchema.addField('ref.'+keyName, type=typeStr)
    for keyName in srcSchema.getNames():
        field = srcSchema.find(keyName).field
        typeStr = field.getTypeString()
        mergedSchema.addField('src.'+keyName, type=typeStr)
                
    mergedCatalog = afwTable.BaseCatalog(mergedSchema)

    refKeys = []
    for keyName in refSchema.getNames():
        refKeys.append((refSchema.find(keyName).key, mergedSchema.find('ref.' + keyName).key))
    srcKeys = []
    for keyName in srcSchema.getNames():
        srcKeys.append((srcSchema.find(keyName).key, mergedSchema.find('src.' + keyName).key))

    for match in matchlist:
        record = mergedCatalog.addNew()
        for key in refKeys:
            keyIn = key[0]
            keyOut = key[1]
            record.set(keyOut, match.first.get(keyIn))
        for key in srcKeys:
            keyIn = key[0]
            keyOut = key[1]
            record.set(keyOut, match.second.get(keyIn))

    if butler is not None:
        butler.put(mergedCatalog, 'matchedList')
        return 
    elif fileName is not None:
        mergedCatalog.writeFits(fileName)
        return
    else:
        print '** No file has been written.'
        return
