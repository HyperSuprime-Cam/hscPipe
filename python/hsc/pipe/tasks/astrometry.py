#!/usr/bin/env python

import lsst.pex.config as pexConfig
import lsst.afw.cameraGeom as afwCG
import lsst.pipe.base as pipeBase
import lsst.meas.astrom as measAstrom
import lsst.pipe.tasks.astrometry as ptAstrometry
import hsc.meas.astrom.astrom as hscAstrom

class HscAstrometryConfig(ptAstrometry.AstrometryConfig):
    solver = pexConfig.ConfigField(
        dtype=hscAstrom.TaburAstrometryConfig,
        doc = "Configuration for the Tabur astrometry solver"
        )


# Use hsc.meas.astrom, failing over to lsst.meas.astrom
class HscAstrometryTask(ptAstrometry.AstrometryTask):
    ConfigClass = HscAstrometryConfig
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


    def distort(self, exposure, sources):
        if exposure.getWcs().hasDistortion():
            for s in sources:
                s.set(self.centroidKey, s.getCentroid())
            return (0,0), (exposure.getWidth(), exposure.getHeight())
        return super(HscAstrometryTask, self).distort(exposure, sources)

