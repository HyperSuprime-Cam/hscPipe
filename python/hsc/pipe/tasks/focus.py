import numpy as np

from lsst.pex.config import Config, Field, ListField

import lsst.afw.math as afwMath
import lsst.afw.geom as afwGeom
import lsst.afw.cameraGeom as afwCG
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


# get solution of the corrrelation function (reconstructed_focus -> true_focus) by Newton's method
def getCorrectedFocus(f, df, corr_coefficients, epsilon=1e-4, n=1000):
    ff = 0.001
    corr = np.poly1d(corr_coefficients)
    dcorr_dtruefocus = corr.deriv()
    i = 0
    while i < n:
        fff = ff - (corr(ff)-f)/dcorr_dtruefocus(ff)
        i += 1
        if np.abs((fff-ff)/ff) < epsilon:
            return fff, np.abs(1./dcorr_dtruefocus(fff))*df
        ff = fff
    raise RuntimeError("Cannot solve for corrected focus: %s %s %s --> %d %f" %
                       (f, df, corr_coefficients, i, (fff-ff/ff)))

def getFocusOffset(ccd, config):
    # physical offsets are +/-0.2 mm, but offsets in the hexapod coordinates are +/-0.12 mm, which is derived from ZEMAX simulations.
    if ccd in config.aboveList:
        return config.offset
    if ccd in config.belowList:
        return -1*config.offset
    raise KeyError("CCD identifier %s not in configuration: %s %s" % (ccd, config.aboveList, config.belowList))

def getDistanceFromFocus(d_icSrc, d_ccdObj, d_dims, zemaxFilename, config, plot_filename=None):
    # set up radial bins. Bins are optimized so that the area of each bin is the same.
    l_binr_e = config.radialBinEdges # [16600, 17380.580580580579, 17728.128128128126, 18000]
    l_binr = config.radialBinCenters # [17112.514461756149, 17563.380665628181, 17868.148132145379]
    l_binr_le = l_binr_e[0:-1]
    l_binr_ue = l_binr_e[1:]

    # make selection on data and get rms^2 for each bin, ccd by ccd
    d_l_rmssq = dict() # rmssq list for radial bin, which is dictionary for each ccd

    for ccd in d_icSrc:
        # use only objects classified as PSF candidate
        icSrc = d_icSrc[ccd][d_icSrc[ccd].get("calib.psf.candidate")]

        # prepare for getting distance from center for each object
        ccdObj = d_ccdObj[ccd]
        x1, y1 = d_dims[ccd]
        # getMm() currently provides a number in pixel.
        # Note that we constructed the zemax values alpha(r), d0(r), and this r is in pixel.
        # We should update this code when getMm() is updated.
        u_llc, v_llc = ccdObj.getPositionFromPixel(afwGeom.PointD(0., 0.)).getMm()
        u_lrc, v_lrc = ccdObj.getPositionFromPixel(afwGeom.PointD(x1, 0.)).getMm()
        u_ulc, v_ulc = ccdObj.getPositionFromPixel(afwGeom.PointD(0., y1)).getMm()
        u_urc, v_urc = ccdObj.getPositionFromPixel(afwGeom.PointD(x1, y1)).getMm()

        l_r = list()
        l_rmssq = list()
        for _icSrc in icSrc:
            # reject blended objects
            if len(_icSrc.getFootprint().getPeaks()) != 1:
#                print "reject a blended object at CCD:%s ,(%f, %f)" % (ccd, _icSrc.getX(), _icSrc.getY())
                continue

            # calculate distance from center for each objects
            x = _icSrc.getX()
            y = _icSrc.getY()

            u_l = (u_lrc-u_llc)/x1*x+u_llc
            u_u = (u_urc-u_ulc)/x1*x+u_ulc
            u = (u_u-u_l)/y1*y+u_l

            v_l = (v_lrc-v_llc)/x1*x+v_llc
            v_u = (v_urc-v_ulc)/x1*x+v_ulc
            v = (v_u-v_l)/y1*y+v_l
            l_r.append(np.sqrt(u**2 + v**2))

            # calculate rms^2
            mom = _icSrc.get(config.shape)
            l_rmssq.append((mom.getIxx() + mom.getIyy())*config.pixelScale**2) # convert from pixel^2 to mm^2

        # calculate median rms^2 for each radial bin
        l_r = np.array(l_r)
        l_rmssq = np.array(l_rmssq)
        l_rmssq_median = list()
        for binr_le, binr_ue in zip(l_binr_le, l_binr_ue):
            sel = np.logical_and(l_r > binr_le, l_r < binr_ue)
            l_rmssq_median.append(np.median(l_rmssq[sel]))
        d_l_rmssq[ccd] = np.ma.masked_array(l_rmssq_median, mask = np.isnan(l_rmssq_median))

    # get ZEMAX values
    d = np.loadtxt(zemaxFilename)

    interpStyle = afwMath.stringToInterpStyle("NATURAL_SPLINE")
    s_alpha = afwMath.makeInterpolate(d[:,0], d[:,1], interpStyle).interpolate
    s_d0 = afwMath.makeInterpolate(d[:,0], d[:,2], interpStyle).interpolate

    # calculate rms^2 for each CCD pair
    l_ccd_pair = zip(config.belowList, config.aboveList)
    l_l_d = list()
    for ccd_pair in l_ccd_pair:
        l_d = list()
        for i_binr, binr in enumerate(l_binr):
            rmssq_p = d_l_rmssq[ccd_pair[1]][i_binr]
            rmssq_m = d_l_rmssq[ccd_pair[0]][i_binr]
            rmssq_diff = rmssq_p - rmssq_m
            delta = getFocusOffset(ccd_pair[1], config)
            alpha = s_alpha(binr)
            d = rmssq_diff/4./alpha/delta + s_d0(binr)
            l_d.append(d)
        l_l_d.append(np.array(l_d))

    l_l_d = np.ma.masked_array(l_l_d, mask = np.isnan(l_l_d))
    d_reconstructed = np.ma.median(l_l_d)
    d_npoint = np.sum(np.invert(l_l_d.mask))
    d_reconstructed_err = np.ma.std(l_l_d)*np.sqrt(np.pi/2.)/np.sqrt(d_npoint)

    if config.doPlot == True:
        if not plot_filename:
            raise ValueError("no filename for focus plot")
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        l_marker = ["o", "x", "d", "^", "<", ">"]
        l_color = ["blue", "green", "red", "cyan", "magenta", "yellow"]
        for i_ccdpair, ccd_pair in enumerate(l_ccd_pair):
            delta_plot = np.ma.masked_array([getFocusOffset(ccd_pair[0], config),
                                             getFocusOffset(ccd_pair[1], config)])
            rmssq_plot = np.ma.masked_array([d_l_rmssq[ccd_pair[0]], d_l_rmssq[ccd_pair[1]]])
            for i_binr in range(len(l_binr)):
                plt.plot(delta_plot, rmssq_plot[:, i_binr], "%s--" % l_marker[i_ccdpair], color = l_color[i_binr])
        plt.savefig(plot_filename)

    d_corrected, d_corrected_err = getCorrectedFocus(d_reconstructed, d_reconstructed_err, config.corrCoeff)
    return d_corrected, d_corrected_err, d_reconstructed, d_reconstructed_err


def run(rerun, frame, doPlot = False):
    l_ccd = [104, 105, 106, 107, 108, 109, 110, 111]

    butler = hscCamera.getButler("HSC", rerun)
    config = FocusConfig()
    config.update(**hscFocusConfig.toDict()) # copying hscConfig
    config.doPlot = doPlot
    config.freeze()

    # get filter, size
    import pyfits
    try:
        header = pyfits.getheader(butler.get('calexp_filename', dataId = {"visit": frame, "ccd": 104})[0])
    except IOError:
        raise
    filter = header["FILTER"]

    # get data and exposure
    d_icSrc = dict()
    d_dims = dict()
    d_ccd = dict()
    for ccd in l_ccd:
        try:
            d_icSrc[ccd] = butler.get('icSrc', dataId = {"visit": int(frame), "ccd": ccd})
        except RuntimeError:
            raise
        exposure = butler.get('calexp', dataId = {"visit": int(frame), "ccd": ccd})
        d_dims[ccd] = exposure.getDimensions()
        d_ccd[ccd] = afwCG.cast_Ccd(exposure.getDetector())

    plot_filename = "focus_%s.png" % frame

    if filter == "g":
        zemax_config = 9
    elif filter == "r":
        zemax_config = 1
    elif filter == "i":
        zemax_config = 3
    elif filter == "z":
        zemax_config = 5
    elif filter == "y":
        zemax_config = 7
    zemax_filename = os.path.join(os.environ["OBS_SUBARU_DIR"], "hsc", "zemax_config%s_0.0.dat" % zemax_config)

    return getDistanceFromFocus(d_icSrc, d_ccd, d_dims, zemax_filename, config,
                                plot_filename=plot_filename)

if __name__ == "__main__":
    import sys
    rerun = "miyatake-focus-test"
    frame = int(sys.argv[1])
    doPlot = True
    d, d_err, d_uncorr, d_uncorr_err = run(rerun, frame, doPlot = doPlot)
    print "focus before correction: %s +/- %s" % (d_uncorr, d_uncorr_err)
    print "focus after correction: %s +/- %s" % (d, d_err)
