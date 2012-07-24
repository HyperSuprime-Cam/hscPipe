
class QaConfig(pexConfig.Config):
    measureSeeing = pexConfig.ConfigField(dtype=hscSeeing.MeasureSeeingConfig, doc="Measure seeing")


class QaTask(pexConfig.Config):
    ConfigClass = QaConfig
    def __init__(self, *args, **kwargs):
    super(QaTask, self).__init(*args, **kwargs)
    self.makeSubtask("measureSeeing", hscSeeing.MeasureSeeingTask)

    def run(self, dataRef, exposure, matches):
        metadata = exposure.getMetadata()

        # = flags
        # this part should be done by calculating merit functions somewhere else in a polite manner.
        metadata.set('FLAG_AUTO', 0)
        metadata.set('FLAG_USR', 0)
        metadata.set('FLAG_TAG', 1)

        # = info for data management
        # = rerun; assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'

        rerunName = self.getRerunName(dataRef)
        metadata.set('RERUN', rerunName)


        # === For QA db
        sdssFluxes = retrieveReferenceMagnitudes(calib.matches)
        #print '*** returned sdssFluxes:', sdssFluxes
        #writeMatchesToBintableFits(calib.matches, butler=sensorRef)
        writeEnhancedMatchesToBintableFits(calib.matches, calib.matchMeta, refMags=sdssFluxes, butler=sensorRef)

    def getRerunName(self, sensorRef):
        return sensorRef.getButler().mapper.rerun




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
