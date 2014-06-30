import numpy as np

from lsst.pex.config import Config, Field, ListField

import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as afwCG
import lsst.pex.exceptions as pexExceptions
import hsc.pipe.base.camera as hscCamera

class FocusConfig(Config):
    # Defaults are appropriate for HSC, but also shouldn't get in the way for Suprime-Cam
    # (because Suprime-Cam CCDs aren't indexed over 10).
    corrCoeff = ListField(dtype=float, default=[8.238421, 1.607829, 1.563773, 0.029580],
                          doc="Correction polynomial coefficients: reconstructed_focus = corr(true_focus)")
    aboveList = ListField(dtype=int, default=[107, 104, 111, 108], doc="Indices of CCDs above focus")
    belowList = ListField(dtype=int, default=[105, 106, 109, 110], doc="Indices of CCDs below focus")
    offset = Field(dtype=float, default=0.12, doc="Focus offset for CCDs")
    radialBinEdges = ListField(dtype=float, default=[16600, 17380.580580580579, 17728.128128128126, 18000],
                               doc="Radii edges for bins")
    radialBinCenters = ListField(dtype=float,
                                 default=[17112.514461756149, 17563.380665628181, 17868.148132145379],
                                 doc="Radii centers for bins")
    doPlot = Field(dtype=bool, default=False, doc="Plot focus calculation?")
    shape = Field(dtype=str, default="shape.sdss", doc="Measurement to use for shape")
    pixelScale = Field(dtype=float, default=0.015, doc="Conversion factor for pixel scale --> mm")

    def isFocusCcd(self, ccdId):
        return ccdId in self.aboveList or ccdId in self.belowList

    def validate(self):
        super(FocusConfig, self).validate()
        numRadialBins = len(self.radialBinCenters)
        if len(self.radialBinEdges) != numRadialBins + 1:
            raise RuntimeError("Expected %d radialBinEdges for the %d radialBinCenters" %
                               (numRadialBins + 1, numRadialBins))
        if len(self.aboveList) != len(self.belowList):
            raise RuntimeError("List of CCDs above and below focus not of equal length: %d %d" %
                               (len(self.aboveList), len(self.belowList)))


# get corrected focus. corrCoeff is polynomials of a function that maps a true focus error to a focus error reconstructed which is obtained from calibration data empirically. This function calculates inverse of the correction function by Newton's method.
def getCorrectedFocusError(f, df, corrCoeff, epsilon=1e-4, n=1000):
    ff = 0.001
    corrFunc = np.poly1d(corrCoeff)
    dCorrFuncDTrueFocus = corrFunc.deriv() # derivative of the correction function
    i = 0
    while i < n:
        fff = ff - (corrFunc(ff)-f)/dCorrFuncDTrueFocus(ff)
        i += 1
        if np.abs((fff-ff)/ff) < epsilon:
            return fff, np.abs(1./dCorrFuncDTrueFocus(fff))*df
        ff = fff
    raise RuntimeError("Cannot solve for corrected focus: %s %s %s --> %d %f" %
                       (f, df, corrCoeff, i, (fff-ff/ff)))

def getFocusCcdOffset(ccd, config):
    # physical offsets are +/-0.2 mm, but offsets in the hexapod coordinates are +/-0.12 mm, which is derived from ZEMAX simulations.
    if ccd in config.aboveList:
        return config.offset
    if ccd in config.belowList:
        return -1*config.offset
    raise KeyError("CCD identifier %s not in configuration: %s %s" % (ccd, config.aboveList, config.belowList))

def getDistanceFromFocus(dIcSrc, dCcd, dCcdDims, zemaxFilename, config, plotFilename=None):
    # Focus error is measured by using rms^2 of stars on focus CCDs.
    # If there is a focus error d, rms^2 can be written as
    # rms^2 = rms_atm^2 + rms_opt_0^2 + alpha*d^2,
    # where rms_atm is from atmosphere and rms_opt if from optics with out any focus error. 
    # On the focus CCDs which have +/-delta offset, the equation becomes
    # rms_+^2 = rms_atm^2 + rms_opt_0^2 + alpha(d+delta)^2
    # rms_-^2 = rms_atm^2 + rms_opt_0^2 + alpha(d-delta)^2
    # Thus, the difference of these rms^2 gives the focus error as
    # d = (rms_+^2 - rms_-^2)/(4 alpha delta)
    # alpha is determined by ZEMAX simulations. It turned out that alpha is a function of distance from the center of FOV r.
    # Also the best focus varies as a function of r. Thus the focus error can be rewritten as
    # d(r) = (rms_+(r)^2 - rms_-(r)^2)/(4 alpha(r) delta) + d0(r)
    # I take a pair of CCDs on the corner, divide the focus CCDs into radian bins, calculate focus error d for each radial bin with alpha and d0 values at this radius, and then take median of these focus errors for all the radian bins and CCD pairs.
    # rms^2 is measured by shape.simple. Although I intend to include minimum measurement bias, there exists still some bias. This is corrected by getCorrectedFocusError() at the end, which is a polynomial function derived by calibration data (well-behaved focus sweeps).

    # set up radial bins
    lRadialBinEdges = config.radialBinEdges
    lRadialBinCenters = config.radialBinCenters
    lRadialBinsLowerEdges = lRadialBinEdges[0:-1]
    lRadialBinsUpperEdges = lRadialBinEdges[1:]

    # make selection on data and get rms^2 for each bin, CCD by CCD
    dlRmssq = dict() # rmssq list for radial bin, which is dictionary for each ccd

    for ccdId in dIcSrc.keys():
        # use only objects classified as PSF candidate
        icSrc = dIcSrc[ccdId][dIcSrc[ccdId].get("calib.psf.candidate")]

        # prepare for getting distance from center for each object
        ccd = dCcd[ccdId]
        x1, y1 = dCcdDims[ccdId]
        # getMm() currently provides a number in pixel.
        # Note that we constructed the zemax values alpha(r), d0(r), and this r is in pixel.
        # We should update this code when getMm() is updated.
        uLlc, vLlc = ccd.getPositionFromPixel(afwGeom.PointD(0., 0.)).getMm()
        uLrc, vLrc = ccd.getPositionFromPixel(afwGeom.PointD(x1, 0.)).getMm()
        uUlc, vUlc = ccd.getPositionFromPixel(afwGeom.PointD(0., y1)).getMm()
        uUrc, vUrc = ccd.getPositionFromPixel(afwGeom.PointD(x1, y1)).getMm()

        lDistanceFromCenter = list()
        lRmssq = list()
        for s in icSrc:
            # reject blended objects
            if len(s.getFootprint().getPeaks()) != 1:
                continue

            # calculate distance from center for each objects
            x = s.getX()
            y = s.getY()

            uL = (uLrc-uLlc)/x1*x+uLlc
            uU = (uUrc-uUlc)/x1*x+uUlc
            u = (uU-uL)/y1*y+uL

            vL = (vLrc-vLlc)/x1*x+vLlc
            vU = (vUrc-vUlc)/x1*x+vUlc
            v = (vU-vL)/y1*y+vL
            lDistanceFromCenter.append(np.sqrt(u**2 + v**2))

            # calculate rms^2
            mom = s.get(config.shape)
            lRmssq.append((mom.getIxx() + mom.getIyy())*config.pixelScale**2) # convert from pixel^2 to mm^2

        # calculate median rms^2 for each radial bin
        lDistanceFromCenter = np.array(lDistanceFromCenter)
        lRmssq = np.array(lRmssq)
        lRmssqMedian = list()
        for radialBinLowerEdge, radialBinUpperEdge in zip(lRadialBinsLowerEdges, lRadialBinsUpperEdges):
            sel = np.logical_and(lDistanceFromCenter > radialBinLowerEdge, lDistanceFromCenter < radialBinUpperEdge)
            lRmssqMedian.append(np.median(lRmssq[sel]))
        dlRmssq[ccdId] = np.ma.masked_array(lRmssqMedian, mask = np.isnan(lRmssqMedian))

    # get ZEMAX values
    d = np.loadtxt(zemaxFilename)

    interpStyle = afwMath.stringToInterpStyle("NATURAL_SPLINE")
    sAlpha = afwMath.makeInterpolate(d[:,0], d[:,1], interpStyle).interpolate
    sD0 = afwMath.makeInterpolate(d[:,0], d[:,2], interpStyle).interpolate

    # calculate rms^2 for each CCD pair
    lCcdPairs = zip(config.belowList, config.aboveList)
    llFocurErrors = list()
    for ccdPair in lCcdPairs:
        lFocusErrors = list()
        for i, radialBinCenter in enumerate(lRadialBinCenters):
            rmssqAbove = dlRmssq[ccdPair[1]][i]
            rmssqBelow = dlRmssq[ccdPair[0]][i]
            rmssqDiff = rmssqAbove - rmssqBelow
            delta = getFocusCcdOffset(ccdPair[1], config)
            alpha = sAlpha(radialBinCenter)
            focusError = rmssqDiff/4./alpha/delta + sD0(radialBinCenter)
            lFocusErrors.append(focusError)
        llFocurErrors.append(np.array(lFocusErrors))

    llFocurErrors = np.ma.masked_array(llFocurErrors, mask = np.isnan(llFocurErrors))
    reconstructedFocusError = np.ma.median(llFocurErrors)
    n = np.sum(np.invert(llFocurErrors.mask))
    reconstructedFocusErrorStd= np.ma.std(llFocurErrors)*np.sqrt(np.pi/2.)/np.sqrt(n)

    if config.doPlot == True:
        if not plotFilename:
            raise ValueError("no filename for focus plot")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        lMarker = ["o", "x", "d", "^", "<", ">"]
        lColor = ["blue", "green", "red", "cyan", "magenta", "yellow"]
        for i, ccdPair in enumerate(lCcdPairs):
            delta_plot = np.ma.masked_array([getFocusCcdOffset(ccdPair[0], config),
                                             getFocusCcdOffset(ccdPair[1], config)])
            rmssq_plot = np.ma.masked_array([dlRmssq[ccdPair[0]], dlRmssq[ccdPair[1]]])
            for j in range(len(lRadialBinCenters)):
                plt.plot(delta_plot, rmssq_plot[:, j], "%s--" % lMarker[i], color = lColor[j])
        plt.savefig(plotFilename)

    correctedFocusError, correctedFocusErrorStd = getCorrectedFocusError(reconstructedFocusError, reconstructedFocusErrorStd, config.corrCoeff)
    return correctedFocusError, correctedFocusErrorStd, reconstructedFocusError, reconstructedFocusErrorStd


def run(rerun, frameId, doPlot = False):
    lCcdId = [104, 105, 106, 107, 108, 109, 110, 111]

    butler = hscCamera.getButler("HSC", rerun)
    config = FocusConfig()
    config.update(**hscFocusConfig.toDict()) # copying hscConfig
    config.doPlot = doPlot
    config.freeze()

    # get filter, size
    try:
        metadata = afwImage.readMetadata(butler.get('calexp_filename', dataId = {"visit": frameId, "ccd": 104})[0])
    except pexExceptions.LsstCppException:
        raise
    filter = metadata.get("FILTER")

    # get data and exposure
    dIcSrc = dict()
    dCcdDims = dict()
    dCcd = dict()
    for ccdId in lCcdId:
        try:
            dIcSrc[ccdId] = butler.get('icSrc', dataId = {"visit": int(frameId), "ccd": ccdId})
        except RuntimeError:
            raise
        exposure = butler.get('calexp', dataId = {"visit": int(frameId), "ccd": ccdId})
        dCcdDims[ccdId] = exposure.getDimensions()
        dCcd[ccdId] = afwCG.cast_Ccd(exposure.getDetector())

    plotFilename = "focus_%s.png" % frameId

    if filter == "g":
        zemaxConfig = 9
    elif filter == "r":
        zemaxConfig = 1
    elif filter == "i":
        zemaxConfig = 3
    elif filter == "z":
        zemaxConfig = 5
    elif filter == "y":
        zemaxConfig = 7
    zemaxFilename = os.path.join(os.environ["OBS_SUBARU_DIR"], "hsc", "zemaxConfig%s_0.0.dat" % zemaxConfig)

    return getDistanceFromFocus(dIcSrc, dCcd, dCcdDims, zemaxFilename, config,
                                plotFilename=plotFilename)

if __name__ == "__main__":
    import sys
    rerun = "miyatake-focus-test"
    frameId = int(sys.argv[1])
    doPlot = True
    correctedFocusError, correctedFocusErrorStd, reconstructedFocusError, reconstructedFocusErrorStd = run(rerun, frameId, doPlot = doPlot)
    print "focus before correction: %s +/- %s" % (reconstructedFocusError, reconstructedFocusErrorStd)
    print "focus after correction: %s +/- %s" % (correctedFocusError, correctedFocusErrorStd)
