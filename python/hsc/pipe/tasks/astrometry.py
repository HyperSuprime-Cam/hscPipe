#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.meas.astrom as measAstrom
import lsst.pipe.tasks.astrometry as ptAstrometry
import hsc.meas.astrom.astrom as hscAstrom
##== FH added for QA output
import lsst.meas.astrom.sip as astromSip


class HscAstrometryConfig(ptAstrometry.AstrometryConfig):
    solver = pexConfig.ConfigField(
        dtype=hscAstrom.TaburAstrometryConfig,
        doc = "Configuration for the Tabur astrometry solver"
        )


# Use hsc.meas.astrom, failing over to lsst.meas.astrom
class HscAstrometryTask(ptAstrometry.AstrometryTask):
    ConfigClass = HscAstrometryConfig

##== FH added for QA output
    @pipeBase.timeMethod
    def run(self, exposure, sources):
        """AstrometryTask an exposure: PSF, astrometry and photometry
        
        @param exposure Exposure to calibrate
        @param sources List of measured sources
        @return a pipeBase.Struct with fields:
        - matches: Astrometric matches
        - matchMeta: Metadata for astrometric matches
        """
        assert exposure is not None, "No exposure provided"
        
        llc, size = self.distort(exposure, sources)
        oldCentroidKey = sources.table.getCentroidKey()
        sources.table.defineCentroid(self.centroidKey, sources.table.getCentroidErrKey(),
                                     sources.table.getCentroidFlagKey())
        ##== FH calls astrometryQa() 
        matches, matchMeta = self.astrometryQa(exposure, sources, llc=llc, size=size)
        sources.table.defineCentroid(oldCentroidKey, sources.table.getCentroidErrKey(),
                                     sources.table.getCentroidFlagKey())

        ##== FH calls undistortQa() 
        self.undistortQa(exposure, sources, matches)
        
        
        return pipeBase.Struct(
            matches = matches,
            matchMeta = matchMeta,
            )

    @pipeBase.timeMethod
    def astrometry(self, exposure, sources, llc=(0,0), size=None):
        """Solve astrometry to produce WCS

        @param exposure Exposure to process
        @param sources Sources
        @param llc Lower left corner (minimum x,y)
        @param size Size of exposure
        @return Star matches, match metadata
        """
        assert exposure, "No exposure provided"

        self.log.log(self.log.INFO, "Solving astrometry")

        try:
            import hsc.meas.astrom as hscAst
        except ImportError:
            hscAst = None

        wcs = exposure.getWcs()
        if wcs is None or hscAst is None:
            self.log.log(self.log.WARN, "Unable to use hsc.meas.astrom; reverting to lsst.meas.astrom")
            return ptAstrometry.AstrometryTask.astrometry(exposure, sources, llc=llc, size=size)

        if size is None:
            size = (exposure.getWidth(), exposure.getHeight())

        try:
            astrometer = hscAstrom.TaburAstrometry(self.config.solver, log=self.log)
            astrom = astrometer.determineWcs(sources, exposure)
            if astrom is None:
                raise RuntimeError("hsc.meas.astrom failed to determine the WCS")
        except Exception, e:
            self.log.log(self.log.WARN, "hsc.meas.astrom failed (%s); trying lsst.meas.astrom" % e)
            astrometer = measAstrom.Astrometry(self.config.solver, log=self.log)
            astrom = astrometer.determineWcs(sources, exposure)
        
        if astrom is None:
            raise RuntimeError("Unable to solve astrometry for %s", exposure.getDetector().getId())

        wcs = astrom.getWcs()
        matches = astrom.getMatches()
        matchMeta = astrom.getMatchMetadata()
        if matches is None or len(matches) == 0:
            raise RuntimeError("No astrometric matches for %s", exposure.getDetector().getId())
        self.log.log(self.log.INFO, "%d astrometric matches for %s" % \
                     (len(matches), exposure.getDetector().getId()))
        exposure.setWcs(wcs)

        # Apply WCS to sources
        for source in sources:
            distorted = source.get(self.centroidKey)
            sky = wcs.pixelToSky(distorted.getX(), distorted.getY())
            source.setCoord(sky) 

        self.display('astrometry', exposure=exposure, sources=sources, matches=matches)

        return matches, matchMeta

##== FH added for QA output
    def astrometryQa(self, exposure, sources, llc=(0,0), size=None):
        """Solve astrometry to produce WCS

        @param exposure Exposure to process
        @param sources Sources
        @param llc Lower left corner (minimum x,y)
        @param size Size of exposure
        @return Star matches, match metadata
        """
        assert exposure, "No exposure provided"

        metadata = exposure.getMetadata()
        
        self.log.log(self.log.INFO, "QA astrometry: Solving astrometry and recording QA outputs")
        
        try:
            import hsc.meas.astrom as hscAst
        except ImportError:
            hscAst = None

        wcs = exposure.getWcs()
        if wcs is None or hscAst is None:
            self.log.log(self.log.WARN, "Unable to use hsc.meas.astrom; reverting to lsst.meas.astrom")
            return ptAstrometry.AstrometryTask.astrometry(exposure, sources, llc=llc, size=size)

        if size is None:
            size = (exposure.getWidth(), exposure.getHeight())

        try:
            astrometer = hscAstrom.TaburAstrometry(self.config.solver, log=self.log)
            astrom = astrometer.determineWcs(sources, exposure)
            if astrom is None:
                raise RuntimeError("hsc.meas.astrom failed to determine the WCS")
        except Exception, e:
            self.log.log(self.log.WARN, "hsc.meas.astrom failed (%s); trying lsst.meas.astrom" % e)
            astrometer = measAstrom.Astrometry(self.config.solver, log=self.log)
            astrom = astrometer.determineWcs(sources, exposure)
        
        if astrom is None:
            raise RuntimeError("Unable to solve astrometry for %s", exposure.getDetector().getId())

        wcs = astrom.getWcs()
        matches = astrom.getMatches()
        matchMeta = astrom.getMatchMetadata()
        if matches is None or len(matches) == 0:
            raise RuntimeError("No astrometric matches for %s", exposure.getDetector().getId())
        self.log.log(self.log.INFO, "%d astrometric matches for %s" % \
                     (len(matches), exposure.getDetector().getId()))
        exposure.setWcs(wcs)

        # Apply WCS to sources
        for source in sources:
            distorted = source.get(self.centroidKey)
            sky = wcs.pixelToSky(distorted.getX(), distorted.getY())
            source.setCoord(sky) 

        self.display('astrometry', exposure=exposure, sources=sources, matches=matches)

        ##== FH added for QA output
        metadata.set('NOBJ_BRIGHT', len(sources))
        metadata.set('NOBJ_MATCHED', len(matches))  # N matched objects ; This should be num of all matched objects in hscAstrom.match()!!
        metadata.set('WCS_NOBJ', len(matches))  # N matched objects used to solve astrometry

        return matches, matchMeta

##== FH added for QA output
    @pipeBase.timeMethod
    def undistortQa(self, exposure, sources, matches):
        """Undistort matches after solving astrometry, resolving WCS

        @param exposure Exposure of interest
        @param sources Sources on image (no distortion applied)
        @param matches Astrometric matches
        @param distortion Distortion model
        """
        assert exposure, "No exposure provided"
        assert sources, "No sources provided"
        assert matches, "No matches provided"

        # Undo distortion in matches
        self.log.log(self.log.INFO, "QA astrometry: Removing distortion correction.")

        metadata = exposure.getMetadata()
        #metadata.set('NOBJ_BRIGHT', len(sources)) # now set in astrometryQa()

        # dummy values in the case of 'without distortion'
        metadata.set('WCS_SIPORDER', 0)
        metadata.set('WCS_RMS', -1)        
                     
        # Re-fit the WCS with the distortion undone
        if self.config.solver.calculateSip:
            self.log.log(self.log.INFO, "QA astrometry: Refitting WCS with distortion removed")
            sip = astromSip.CreateWcsWithSip(matches, exposure.getWcs(), self.config.solver.sipOrder)
            wcs = sip.getNewWcs()
            self.log.log(self.log.INFO, "QA astrometry: Astrometric scatter: %f arcsec (%s non-linear terms)" %
                         (sip.getScatterOnSky().asArcseconds(), "with" if wcs.hasDistortion() else "without"))
            exposure.setWcs(wcs)

            # metadata for QA
            if wcs.hasDistortion() is True:
                sipOrder = self.config.solver.sipOrder
            else: sipOrder = 0
            metadata.set('WCS_SIPORDER', self.config.solver.sipOrder)
            metadata.set('WCS_RMS', sip.getScatterOnSky().asArcseconds()) # arcsec rms
            
            # Apply WCS to sources
            for index, source in enumerate(sources):
                sky = wcs.pixelToSky(source.getX(), source.getY())
                source.setCoord(sky)
        else:
            self.log.log(self.log.WARN, "Not calculating a SIP solution; matches may be suspect")
        
        self.display('astrometry', exposure=exposure, sources=sources, matches=matches)

