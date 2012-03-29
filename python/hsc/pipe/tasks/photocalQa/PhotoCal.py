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
import math, os, sys
import numpy as np

import lsst.meas.astrom as measAst
import lsst.meas.astrom.sip as sip
import lsst.meas.algorithms.utils as malgUtil
import lsst.pex.logging as pexLog
import lsst.pex.config as pexConf
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from lsst.meas.photocal.PhotometricMagnitude import PhotometricMagnitude

try:
#    import matplotlib 
#    matplotlib.use('Agg')
    import matplotlib.pyplot as pyplot
except ImportError:
    pyplot = None

class PhotoCalConfig(pexConf.Config):

    magLimit = pexConf.Field(dtype=float, doc="Don't use objects fainter than this magnitude", default=22.0)
    outputField = pexConf.Field(
        dtype=str, optional=True, default="classification.photometric",
        doc="Name of the flag field that is set for sources used in photometric calibration"
        )
    fluxField = pexConf.Field(
        dtype=str, default="flux.psf", optional=False,
        doc="Name of the source flux field to use.  The associated flag field\n"\
            "('<name>.flags') will be implicitly included in badFlags.\n"
        )
    goodFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=[], 
        doc="List of source flag fields that must be set for a source to be used."
        )
    badFlags = pexConf.ListField(
        dtype=str, optional=False,
        default=["flags.pixel.edge", "flags.pixel.interpolated.any", "flags.pixel.saturated.any"], 
        doc="List of source flag fields that will cause a source to be rejected when they are set."
        )

class PhotoCalTask(pipeBase.Task):
    """Calculate the zero point of an exposure given a ReferenceMatchVector.
    """
    ConfigClass = PhotoCalConfig

    def __init__(self, schema, **kwds):
        """Create the task, pulling input keys from the schema and adding a flag field
        'classification.photometric' that will be set for sources used to determine the
        photometric calibration.
        """
        pipeBase.Task.__init__(self, **kwds)
        self.flux = schema.find(self.config.fluxField).key
        self.fluxErr = schema.find(self.config.fluxField + ".err").key
        self.goodFlags = [schema.find(name).key for name in self.config.goodFlags]
        self.badFlags = [schema.find(self.config.fluxField + ".flags").key]
        self.badFlags.extend(schema.find(name).key for name in self.config.badFlags)
        if self.config.outputField is not None:
            self.output = schema.addField(self.config.outputField, type="Flag",
                                          doc="set if source was used in photometric calibration")
        else:
            self.output = None

    def checkSourceFlags(self, source):
        """Return True if the given source has all good flags set and none of the bad flags set."""
        for k in self.goodFlags:
            if not source.get(k): return False
        for k in self.badFlags:
            if source.get(k): return False
        return True

    @pipeBase.timeMethod
    def selectMatches(self, matches):
        """Select reference/source matches according the criteria specified in the config.
        
        The return value is a ReferenceMatchVector that contains only the selected matches.
        If a schema was passed during task construction, a flag field will be set on sources 
        in the selected matches.

        An exception will be raised if there are no valid matches.

        @param[in] Input ReferenceMatchVector (not modified)
        @return Output ReferenceMatchVector
        """

        self.log.logdebug("Number of input matches: %d" % (len(matches)))
        if len(matches) == 0:
            raise ValueError("No input matches")

        # Only use stars for which the flags indicate the photometry is good.
        afterFlagCut = [m for m in matches if self.checkSourceFlags(m.second)]
        self.log.logdebug("Number of matches after source flag cuts: %d" % (len(afterFlagCut)))
        if len(matches) == 0:
            raise ValueError("All matches eliminated by to source flags")

        refSchema = matches[0].first.schema
        try:
            refKey = refSchema.find("photometric").key
        except:
            self.log.warn("No 'photometric' flag key found in reference schema.")
            refKey = None
        if refKey is not None:
            afterRefCut = [m for m in afterFlagCut if m.first.get(refKey)]
        else:
            afterRefCut = afterFlagCut

        self.log.logdebug("Number of matches after reference catalog cuts: %d" % (len(afterRefCut)))
        if len(afterRefCut) == 0:
            raise RuntimeError("No sources remaining in match list after reference catalog cuts.")
        
        fluxKey = refSchema.find("flux").key
        if self.config.magLimit is not None:
            fluxLimit = 10.0**(-self.config.magLimit / 2.5)
            afterMagCut = [m for m in afterRefCut if (m.first.get(fluxKey) > fluxLimit
                                                      and m.second.get(self.flux) > 0.0)]
        else:
            afterMagCut = [m for m in afterRefCut if m.second.get(self.flux) > 0.0]
        self.log.logdebug("Number of matches after magnitude limit cuts: %d" % (len(afterMagCut)))
        if len(afterRefCut) == 0:
            raise RuntimeError("No sources remaining in match list after magnitude limit cuts.")

        result = afwTable.ReferenceMatchVector()
        for m in afterMagCut:
            if self.output is not None:
                m.second.set(self.output, True)
            result.append(m)
        return result

    @pipeBase.timeMethod
    def extractMagArrays(self, matches):
        """Extract magnitude and magnitude error arrays from the given matches.

        @param[in]  ReferenceMatchVector object containing reference/source matches
        
        @return Struct containing srcMag, refMag, srcMagErr, refMagErr, and errMag arrays.
        """
        srcFlux = np.array([m.second.get(self.flux) for m in matches])
        srcFluxErr = np.array([m.second.get(self.fluxErr) for m in matches])
        if not np.all(np.isfinite(srcFluxErr)):
            self.log.warn("Source catalog does not have flux uncertainties; using sqrt(flux).")
            srcFluxErr = np.sqrt(srcFlux)
        
        refSchema = matches[0].first.schema
        refFluxKey = refSchema.find("flux").key
        refFlux    = np.array([m.first.get(refFluxKey) for m in matches])
        try:
            refFluxErrKey = refSchema.find("flux.err").key
            refFluxErr = np.array([m.first.get(refFluxErrKey) for m in matches])
        except:
            # Catalogue may not have flux uncertainties; HACK
            self.log.warn("Reference catalog does not have flux uncertainties; using sqrt(flux).")
            refFluxErr = np.sqrt(refFlux)

        srcMag = -2.5 * np.log10(srcFlux)
        refMag = -2.5 * np.log10(refFlux)

        # Fitting with error bars in both axes is hard, so transfer all
        # the error to src, then convert to magnitude
        fluxErr = np.hypot(srcFluxErr, refFluxErr)
        magErr = fluxErr / (srcFlux * np.log(10))

        srcMagErr = srcFluxErr / (srcFlux * np.log(10))
        refMagErr = refFluxErr / (refFlux * np.log(10))

        return pipeBase.Struct(
            srcMag = srcMag,
            refMag = refMag,
            magErr = magErr,
            srcMagErr = srcMagErr,
            refMagErr = refMagErr
            )

    @pipeBase.timeMethod
    def run(self, matches):
        """Do photometric calibration - select matches to use and (possibly iteratively) compute
        the zero point.

        @param[in]  matches   Input ReferenceMatchVector (will not be modified).
        
        @return Struct of:
           photocal ---- PhotometricMagnitude object containing the zero point
           arrays ------ Magnitude arrays returned be extractMagArrays
           matches ----- Final ReferenceMatchVector, as returned by selectMathces.
        """
        

        global display, fig
        import lsstDebug
        display = lsstDebug.Info(__name__).display
        if display:
            from matplotlib import pyplot
            try:
                fig.clf()
            except:
                fig = pyplot.figure()

        matches = self.selectMatches(matches)

        arrays = self.extractMagArrays(matches)

        # Fit for zeropoint.  We can run the code more than once, so as to
        # give good stars that got clipped by a bad first guess a second
        # chance.
        # FIXME: these should be config values
        sigma_max = [0.25]                  # maximum sigma to use when clipping
        nsigma = [3.0]                      # clip at nsigma
        useMedian = [True]
        niter = [20]                        # number of iterations

        zp = None                           # initial guess
        for i in range(len(nsigma)):
            r = self.getZeroPoint(arrays.srcMag, arrays.refMag, arrays.magErr, zp0=zp,
                                  useMedian=useMedian[i], sigma_max=sigma_max[i],
                                  nsigma=nsigma[i], niter=niter[i])

            zp = r.zp
            self.log.info("Magnitude zero point: %f +/- %f from %d stars" % (r.zp, r.sigma, r.ngood))

        photocal = PhotometricMagnitude(zeroFlux=1.0, zeroMag=zp)

        self.log.info("QA photocal: Magnitude zero point: %f +/- %f from %d stars" % (r.zp, r.sigma, r.ngood))
        
        ##== FH modified for QA output
        return pipeBase.Struct(
            photocal = photocal,
            arrays = arrays,
            matches = matches,
            zp = r.zp, 
            sigma = r.sigma,
            ngood = r.ngood
            )

    def getZeroPoint(self, src, ref, srcErr=None, zp0=None, 
                     useMedian=True, sigma_max=None, nsigma=2, niter=3):
        """Flux calibration code, returning (ZeroPoint, Distribution Width, Number of stars)

        We perform niter iterations of a simple sigma-clipping algorithm with a a couple of twists:
        1.  We use the median/interquartile range to estimate the position to clip around, and the
        "sigma" to use.
        2.  We never allow sigma to go _above_ a critical value sigma_max --- if we do, a sufficiently
        large estimate will prevent the clipping from ever taking effect.
        3.  Rather than start with the median we start with a crude mode.  This means that a set of magnitude
        residuals with a tight core and asymmetrical outliers will start in the core.  We use the width of
        this core to set our maximum sigma (see 2.)  
        """

        dmag = ref - src
        i = np.argsort(dmag)
        dmag = dmag[i]
        if srcErr is not None:
            dmagErr = srcErr[i]
        else:
            dmagErr = np.ones(len(dmag))

        IQ_TO_STDEV = 0.741301109252802;    # 1 sigma in units of interquartile (assume Gaussian)

        npt = len(dmag)
        ngood = npt
        for i in range(niter):
            if i > 0:
                npt = sum(good)

            center = None
            if i == 0:
                #
                # Start by finding the mode
                #
                nhist = 20
                try:
                    hist, edges = np.histogram(dmag, nhist, new=True)
                except TypeError:
                    hist, edges = np.histogram(dmag, nhist) # they removed new=True around numpy 1.5
                imode = np.arange(nhist)[np.where(hist == hist.max())]

                if imode[-1] - imode[0] + 1 == len(imode): # Multiple modes, but all contiguous
                    if zp0:
                        center = zp0
                    else:
                        center = 0.5*(edges[imode[0]] + edges[imode[-1] + 1])

                    peak = sum(hist[imode])/len(imode) # peak height

                    # Estimate FWHM of mode
                    j = imode[0]
                    while j >= 0 and hist[j] > 0.5*peak:
                        j -= 1
                    j = max(j, 0)
                    q1 = dmag[sum(hist[range(j)])]

                    j = imode[-1]
                    while j < nhist and hist[j] > 0.5*peak:
                        j += 1
                    j = min(j, nhist - 1)
                    j = min(sum(hist[range(j)]), npt - 1)
                    q3 = dmag[j]

                    if q1 == q3:
                        q1 = dmag[int(0.25*npt)]
                        q3 = dmag[int(0.75*npt)]

                    sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

                    if sigma_max is None:
                        sigma_max = 2*sig   # upper bound on st. dev. for clipping. multiplier is a heuristic

                    self.log.logdebug("Photo calibration histogram: center = %.2f, sig = %.2f" 
                                      % (center, sig))

                else:
                    if sigma_max is None:
                        sigma_max = dmag[-1] - dmag[0]

                    center = np.median(dmag)
                    q1 = dmag[int(0.25*npt)]
                    q3 = dmag[int(0.75*npt)]
                    sig = (q3 - q1)/2.3 # estimate of standard deviation (based on FWHM; 2.358 for Gaussian)

            if center is None:              # usually equivalent to (i > 0)
                gdmag = dmag[good]
                if useMedian:
                    center = np.median(gdmag)
                else:
                    gdmagErr = dmagErr[good]
                    center = np.average(gdmag, weights=gdmagErr)

                q3 = gdmag[min(int(0.75*npt + 0.5), npt - 1)]
                q1 = gdmag[min(int(0.25*npt + 0.5), npt - 1)]

                sig = IQ_TO_STDEV*(q3 - q1)     # estimate of standard deviation

            good = abs(dmag - center) < nsigma*min(sig, sigma_max) # don't clip too softly

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            if display:
                from matplotlib import pyplot
                try:
                    fig.clf()

                    axes = fig.add_axes((0.1, 0.1, 0.85, 0.80));

                    axes.plot(ref[good], dmag[good] - center, "b+")
                    axes.errorbar(ref[good], dmag[good] - center, yerr=dmagErr[good], 
                                  linestyle='', color='b')

                    bad = np.logical_not(good)
                    if len(ref[bad]) > 0:
                        axes.plot(ref[bad], dmag[bad] - center, "r+")
                        axes.errorbar(ref[bad], dmag[bad] - center, yerr=dmagErr[bad], 
                                      linestyle='', color='r')

                    axes.plot((-100, 100), (0, 0), "g-")
                    for x in (-1, 1):
                        axes.plot((-100, 100), x*0.05*np.ones(2), "g--")

                    axes.set_ylim(-1.1, 1.1)
                    axes.set_xlim(24, 13)
                    axes.set_xlabel("Reference")
                    axes.set_ylabel("Reference - Instrumental")

                    fig.show()

                    if display > 1:
                        while i == 0 or reply != "c":
                            try:
                                reply = raw_input("Next iteration? [ynhpc] ")
                            except EOFError:
                                reply = "n"

                            if reply == "h":
                                print >> sys.stderr, "Options: c[ontinue] h[elp] n[o] p[db] y[es]"
                                continue

                            if reply in ("", "c", "n", "p", "y"):
                                break
                            else:
                                print >> sys.stderr, "Unrecognised response: %s" % reply

                        if reply == "n":
                            break
                        elif reply == "p":
                            import pdb; pdb.set_trace()
                except Exception, e:
                    print >> sys.stderr, "Error plotting in PhotoCal.getZeroPoint: %s" % e

            #-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

            old_ngood = ngood
            ngood = sum(good)
            if ngood == 0:
                msg = "PhotoCal.getZeroPoint: no good stars remain"

                if i == 0:                  # failed the first time round -- probably all fell in one bin
                    center = np.average(dmag, weights=dmagErr)
                    msg += " on first iteration; using average of all calibration stars"


                self.log.log(self.log.WARN, msg)

                return center, sig, len(dmag)
            elif ngood == old_ngood:
                break

            if False:
                ref = ref[good]
                dmag = dmag[good]
                dmagErr = dmagErr[good]

        dmag = dmag[good]
        dmagErr = dmagErr[good]
        return pipeBase.Struct(
            zp = np.average(dmag, weights=dmagErr),
            sigma = np.std(dmag, ddof=1),
            ngood = len(dmag)
            )
