import os
import sys
import numpy as np
import pyfits
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()

from lsst.pex.config import Config, Field, ListField
import lsst.daf.base as dafBase
import lsst.afw.coord as afwCoord
import lsst.daf.persistence as dafPersist
from lsst.obs.suprimecam import SuprimecamMapper as Mapper
import hsc.pipe.base.camera as hscCamera
import lsst.afw.display.ds9 as ds9
import lsst.afw.cameraGeom as afwCG
import lsst.afw.geom as afwGeom
import focus as focusCode
import lsst.afw.math as afwMath
import lsst.pex.exceptions as pexExceptions
import lsst.afw.image as afwImage

class FocusSweepConfig(Config):
    # Defaults are appropriate for HSC, but also shouldn't get in the way for Suprime-Cam
    # (because Suprime-Cam CCDs aren't indexed over 10).
    corrCoeff = ListField(dtype=float, default=[8.238421, 1.607829, 1.563773, 0.029580],
                          doc="Correction polynomial coefficients: reconstructed_focus = corr(true_focus)")
    aboveList = ListField(dtype=int, default=[107, 104, 111, 108], doc="Indices of CCDs above focus")
    belowList = ListField(dtype=int, default=[105, 106, 109, 110], doc="Indices of CCDs below focus")
    sweepPYList = ListField(dtype=int, default=[105, 107, 108, 110], doc="Indices of CCDs where y becomes larger when sweeping")
    sweepMYList = ListField(dtype=int, default=[104, 106, 110, 111], doc="Indices of CCDs where y becomes smaller when sweeping")
    centerCcdId = Field(dtype=int, default=49, doc="Ceter chip to determine focus from focus sweep")
    offset = Field(dtype=float, default=0.12, doc="Focus offset for CCDs")
    tolXFocusChip = Field(dtype=float, default=15., doc="Tolerance for selecting an object in a focus sweep in x-derection in focus chips [pixel]")
    tolYFocusChip = Field(dtype=float, default=15., doc="Tolerance for selecting an object in a focus sweep in y-derection in focus chips [pixel]")
    tolXCenterChip = Field(dtype=float, default=3., doc="Tolerance for selecting an object in a focus sweep in x-derection in center chips [pixel]")
    tolYCenterChip = Field(dtype=float, default=3., doc="Tolerance for selecting an object in a focus sweep in y-derection in center chips [pixel]")
    dxFocusSweepCenterChip = Field(dtype=float, default=30., doc="Interval of focus sweep in x at center chips[pixel]")
    dxFocusSweepFocusChip = Field(dtype=float, default=34., doc="Interval of focus sweep in x at focus Chips[pixel]")
    dyFocusSweepFocusChip = Field(dtype=float, default=2., doc="Interval of focus sweep in |y| at focus Chips [pixel]")
    radialBinEdges = ListField(dtype=float, default=[16600, 17380.580580580579, 17728.128128128126, 18000],
                               doc="Radii edges for bins")
    radialBinCenters = ListField(dtype=float,
                                 default=[17112.514461756149, 17563.380665628181, 17868.148132145379],
                                 doc="Radii centers for bins")
    doPlot = Field(dtype=bool, default=False, doc="Plot focus calculation?")
    shape = Field(dtype=str, default="shape.simple", doc="Measurement to use for shape")
    pixelScale = Field(dtype=float, default=0.015, doc="Conversion factor for pixel scale --> mm")

    def isFocusOrCenterCcd(self, ccdId):
        return ccdId in self.aboveList or ccdId in self.belowList or ccdId == self.centerCcdId

    def validate(self):
        super(FocusSweepConfig, self).validate()
        numRadialBins = len(self.radialBinCenters)
        if len(self.radialBinEdges) != numRadialBins + 1:
            raise RuntimeError("Expected %d radialBinEdges for the %d radialBinCenters" %
                               (numRadialBins + 1, numRadialBins))
        if len(self.aboveList) != len(self.belowList):
            raise RuntimeError("List of CCDs above and below focus not of equal length: %d %d" %
                               (len(self.aboveList), len(self.belowList)))

def getFocusCcdOffset(ccd, config):
    # physical offsets are +/-0.2 mm, but offsets in the hexapod coordinates are +/-0.12 mm, which is derived from ZEMAX simulations.
    if ccd in config.aboveList:
        return config.offset
    if ccd in config.belowList:
        return -1*config.offset
    raise KeyError("CCD identifier %s not in configuration: %s %s" % (ccd, config.aboveList, config.belowList))


def getFocusSweepIds(src, tolX, tolY, dx = 30., dy = 0.):
    # Find focus sweep using sources do not have children.
    # After selecting focus sweep candidates, we require at least one object in the foscus sweep is PSF candidate.
    #selNoChild = np.ones(src.get("id").shape, dtype = np.bool) # (src.get("deblend.nchild") == 0)
    selNoChild = (src.get("deblend.nchild") == 0)
    srcId = src.get("id")[selNoChild]
    srcX = src.getX()[selNoChild]
    srcY = src.getY()[selNoChild]
    srcPsfCand = src.get("calib.psf.candidate")[selNoChild]
    lFocusSweepId = list()
    slope = dy/dx
    for s in src[selNoChild]:
#        selY = np.abs(s.getY() - srcY) < tolY
        selY = np.abs(s.getY() + slope * s.getX()- srcY) < tolY
        if selY.sum() < 7:
            continue
        flag = False
        lId = list()
        for i in range(7):
            if i == 0:
                criteria = np.abs(s.getX() - srcX[selY]) < tolX
            elif i == 1:
                criteria = np.abs(prev + 2.*dx - srcX[selY]) < tolX
            else:
                criteria = np.abs(prev + dx - srcX[selY]) < tolX
            if criteria.sum() == 1:
                arg = int(np.where(criteria == True)[0][0])
                prev = srcX[selY][arg]
                lId.append(srcId[selY][arg])
            else:
                flag = True
                break
        nPsfCand = 0
        for id in lId:
            arg = np.where(srcId == id)[0]
            nPsfCand += srcPsfCand[arg][0]
        if nPsfCand == 0:
            flag = True
        if flag:
            continue
        lFocusSweepId.append(lId)
    return lFocusSweepId


def determineFocusFromFocusSweep(src, lZ, config, returnData = False):
    print "# of src center CCD: %s %s" % (len(src), np.sum(src.get("calib.psf.candidate")))
    tolX = config.tolXCenterChip # tolerance to pick up focus sweep [pixel]
    tolY = config.tolYCenterChip # tolerance to pick up focus sweep [pixel]
    dx = config.dxFocusSweepCenterChip # dx of focus sweep [pixel]

    lFocusSweepId = getFocusSweepIds(src, tolX, tolY, dx = dx)

    # get rmssq
    lRmssqAll = list()
    for focusSweepId in lFocusSweepId:
        lRmssq = list()
        for id in focusSweepId:
            arg = np.where(src.get("id") == id)[0][0]
            mom = src[arg].get(config.shape)
            lRmssq.append((mom.getIxx() + mom.getIyy())*config.pixelScale**2)
        lRmssqAll.append(lRmssq)

    lRmssqMedian = np.median(lRmssqAll, axis = 0)
    p = np.polyfit(lZ, lRmssqMedian, deg = 2)
    a = p[0]; b = p[1]; c = p[2]
    min = -(b**2-4.*a*c)/4./a
    d0 = -b/2./a
    poly = np.poly1d(p)

    if returnData:
        return d0, lRmssqMedian, lRmssqAll, poly
    else:
        return d0

def testFocusDeterminationByFocusSweep(dSrc, dCcd, dCcdDims, lZ, zemaxFilename, frameId,config, plotFilename=None, debug = False):

    centerCcdId = config.centerCcdId
    lCcdId = ([ccdId for ccdId in config.aboveList] + [ccdId for ccdId in config.belowList])
    lCcdIdPair = zip(config.belowList, config.aboveList)

    # setup bins
    lRadialBinEdges = config.radialBinEdges
    lRadialBinCenters = config.radialBinCenters
    lRadialBinLowerEdges = lRadialBinEdges[0:-1]
    lRadialBinUpperEdges = lRadialBinEdges[1:]

    focusZ, lRmssqMedian, lRmssqAll, poly = determineFocusFromFocusSweep(dSrc[centerCcdId], lZ, config, returnData = True)

    # plot focus sweep
    plt.figure(0, figsize = (16, 8))
    plt.subplot(1,2,1)
    for lRmssq in lRmssqAll:
        plt.plot(lZ, lRmssq, "-", alpha = 0.3)
    plt.plot(lZ, lRmssqMedian, "k--", linewidth = 2)
    lZFine = np.linspace(lZ[0], lZ[-1])
    plt.plot(lZFine, poly(lZFine), "-k", linewidth = 2)
    plt.axvline(focusZ, linestyle = "--", color = "red")
    plt.ylim(0., 0.01)
    plt.xlabel("z[mm]")
    plt.ylabel(r"rms$^2$[mm$^2$]")

    tolX = config.tolXFocusChip # tolerance to pick up focus sweep [pixel]
    tolY = config.tolYFocusChip # tolerance to pick up focus sweep [pixel]

    dlRmssq = dict()
    iCcd = 0
    for ccdId in lCcdId:
        dlRmssq[ccdId] = list()
        src = dSrc[ccdId]

        if debug:
            plt.figure(1)
            plt.subplot(2,4,iCcd+1)
            iCcd += 1
            plt.plot(src.getX(), src.getY(), ".")
            plt.plot(src.getX()[src.get("calib.psf.candidate") == True], src.getY()[src.get("calib.psf.candidate") == True], ".r")

        print "# of src CCD %s: %s %s" % (ccdId, len(src), np.sum(src.get("calib.psf.candidate")))
        try:
            sel = src.get("calib.psf.candidate")
        except RuntimeError, err:
            print "Error", err
            for i in range(len(lRadialBinCenters)):
                dlRmssq[ccdId].append(np.ma.masked_array(0., mask = True))
            continue

        dx = config.dxFocusSweepFocusChip # dx of focus sweep [pixel]
        dy = config.dyFocusSweepFocusChip # dy of focus sweep [pixel]
        if ccdId in config.sweepMYList:
            dy *= - 1

        lFocusSweepIds = getFocusSweepIds(src, tolX, tolY, dx = dx)

        # prepare for getting r
        ccd = dCcd[ccdId]
        x1, y1 = dCcdDims[ccdId]
        u_llc, v_llc = ccd.getPositionFromPixel(afwGeom.PointD(0., 0.)).getMm()
        u_lrc, v_lrc = ccd.getPositionFromPixel(afwGeom.PointD(x1, 0.)).getMm()
        u_ulc, v_ulc = ccd.getPositionFromPixel(afwGeom.PointD(0., y1)).getMm()
        u_urc, v_urc = ccd.getPositionFromPixel(afwGeom.PointD(x1, y1)).getMm()

        # get rmssq and distance from center for each object
        lRmssqAll = list()
        lRAll = list()
        for focusSweepIds in lFocusSweepIds:
            lRmssq = list()
            lR = list()
            lX = list()
            lY = list()
            for id in focusSweepIds:
                arg = np.where(src.get("id") == id)[0][0]

                mom = src[arg].get(config.shape)
                lRmssq.append((mom.getIxx() + mom.getIyy())*config.pixelScale**2)

                x = src[arg].getX()
                y = src[arg].getY()
                lX.append(x)
                lY.append(y)
                u_l = (u_lrc-u_llc)/x1*x+u_llc
                u_u = (u_urc-u_ulc)/x1*x+u_ulc
                u = (u_u-u_l)/y1*y+u_l
                v_l = (v_lrc-v_llc)/x1*x+v_llc
                v_u = (v_urc-v_ulc)/x1*x+v_ulc
                v = (v_u-v_l)/y1*y+v_l
                lR.append(np.sqrt(u**2 + v**2))
            lRmssqAll.append(lRmssq)
            lRAll.append(lR)
            if debug:
                plt.plot(lX, lY, "xk")
                plt.xlim(0, 2048)
                plt.ylim(0, 4096)
                plt.title("%s" % ccdId)

        lRAll = np.array(lRAll)

        for radialBinLowerEdge, radialBinUpperEdge in zip(lRadialBinLowerEdges, lRadialBinUpperEdges):
            selR = np.logical_and(lRAll > radialBinLowerEdge, lRAll < radialBinUpperEdge)
            dlRmssq[ccdId].append(np.ma.median(np.ma.masked_array(lRmssqAll, mask = np.invert(selR)), axis = 0))

    d = np.loadtxt(zemaxFilename)
    interpStyle = afwMath.stringToInterpStyle("NATURAL_SPLINE")
    fa = afwMath.makeInterpolate(d[:,0], d[:,1], interpStyle).interpolate
    fd0 = afwMath.makeInterpolate(d[:,0], d[:,2], interpStyle).interpolate

    lDAll = list()
    for ccdPair in lCcdIdPair:
        if len(dlRmssq[ccdPair[1]][0].shape) == 0 or len(dlRmssq[ccdPair[0]][0].shape) == 0:
            print "pair (%s, %s) is skipped because of no stars in either chip." % (ccdPair[0], ccdPair[1])
            continue
        for i, binr in enumerate(lRadialBinCenters):
            lD = list()
            for j in range(len(lZ)):
                rmssq_p = dlRmssq[ccdPair[1]][i][j]
                rmssq_m = dlRmssq[ccdPair[0]][i][j]
                rmssq_diff = rmssq_p - rmssq_m
                delta = getFocusCcdOffset(ccdPair[1], config)
                a = fa(binr)
                lD.append(rmssq_diff/4./a/delta)
            lD = lD + focusZ + fd0(binr)
            lD = np.array(lD)
            lDAll.append(lD)
    lDAll = np.array(lDAll)
    lDAll = np.ma.masked_array(lDAll, mask = np.isnan(lDAll))
    lDMedian = np.ma.median(lDAll, axis = 0)
    lDNpoint = np.sum(np.invert(lDAll.mask), axis = 0)
    lDStd = np.ma.std(lDAll, axis = 0)*np.sqrt(np.pi/2.)/np.sqrt(lDNpoint)
    lCorrDMedian = list()
    lCorrDStd = list()
    for dMedian, dStd in zip(lDMedian, lDStd):
        corrDMedian, corrDStd = focusCode.getCorrectedFocusError(dMedian - focusZ, dStd, config.corrCoeff)
        lCorrDMedian.append(corrDMedian + focusZ)
        lCorrDStd.append(corrDStd)
    lCorrDMedian = np.array(lCorrDMedian)
    lCorrDStd = np.array(lCorrDStd)

    # fit uncorrected focus error
    S = np.sum(lDStd**-2)
    Sx = np.sum(lZ*lDStd**-2)
    Sy = np.sum(lDMedian*lDStd**-2)
    Sxx = np.sum(lZ**2*lDStd**-2)
    Sxy = np.sum(lZ*lDMedian*lDStd**-2)
    delta = S*Sxx - Sx**2
    a = (Sxx*Sy-Sx*Sxy)/delta
    b = (S*Sxy-Sx*Sy)/delta
    a_err = np.sqrt(Sxx/delta)
    b_err = np.sqrt(S/delta)

    # fit corrected focus error
    S = np.sum(lCorrDStd**-2)
    Sx = np.sum(lZ*lCorrDStd**-2)
    Sy = np.sum(lCorrDMedian*lCorrDStd**-2)
    Sxx = np.sum(lZ**2*lCorrDStd**-2)
    Sxy = np.sum(lZ*lCorrDMedian*lCorrDStd**-2)
    delta = S*Sxx - Sx**2
    a_corr = (Sxx*Sy-Sx*Sxy)/delta
    b_corr = (S*Sxy-Sx*Sy)/delta
    a_corr_err = np.sqrt(Sxx/delta)
    b_corr_err = np.sqrt(S/delta)

    # plot reconstructed z vs input z
    plt.figure(0)
    plt.subplot(1, 2, 2)
    i = 0
    for ccdPair in lCcdIdPair:
        if len(dlRmssq[ccdPair[1]][0].shape) == 0 or len(dlRmssq[ccdPair[0]][0].shape) == 0:
            print "pair (%s, %s) is skipped because of no stars in either chip." % (ccdPair[0], ccdPair[1])
            continue
        for binr in lRadialBinCenters:
            plt.plot(lZ, lDAll[i], "o-", label = "%s, %s" % (ccdPair, binr), alpha = 0.5)
            i += 1
    plt.plot(lZ, b*lZ+a, "k--", linewidth = 2)
    plt.plot(lZ, b_corr*lZ+a_corr, "k--", linewidth = 2)
    plt.suptitle("%s " % frameId)
    plt.errorbar(lZ, lDMedian, lDStd, fmt = "^k", markersize = 10, label = "uncorrected")
    plt.errorbar(lZ, lCorrDMedian, lCorrDStd, fmt = "*k", markersize = 10, label = "corrected")
    plt.plot(lZ, lZ, "r--", linewidth = 1, label = "y=x")
    plt.axhline(focusZ, linestyle = "--", color = "red")
    plt.axvline(focusZ, linestyle = "--", color = "red")
    plt.axis("equal")
    plt.xlabel("z [mm]")
    plt.ylabel("reconstructed z [mm]")
    plt.xlim(lZ[0] - 0.15, lZ[-1] + 0.15)
    plt.ylim(lZ[0] - 0.15, lZ[-1] + 0.15)
    plt.legend(loc = "lower right",  prop = {'size': 10})
    if os.path.dirname(plotFilename) != "" and os.path.exists(os.path.dirname(plotFilename)) == False:
        os.makedirs(os.path.dirname(plotFilename))
    plt.savefig(plotFilename)

def run(rerun, frameId):

    config = FocusSweepConfig()
    config.freeze()

    lCcdId = ([ccdId for ccdId in config.aboveList] + [ccdId for ccdId in config.belowList])
    lCcdId.append(config.centerCcdId)

    butler = hscCamera.getButler("HSC", rerun)


    # get filter, size
    try:
        metadata = afwImage.readMetadata(butler.get('calexp_filename', dataId = {"visit": frameId, "ccd": 104})[0])
    except pexExceptions.LsstCppException:
        raise
    filter = metadata.get("FILTER")
    focusHeader = metadata.get("FOC-VAL")

    # get data and exposure
    dSrc = dict()
    dCcdDims = dict()
    dCcd = dict()
    for ccdId in lCcdId:
        try:
            dSrc[ccdId] = butler.get('src', dataId = {"visit": int(frameId), "ccd": ccdId})
        except RuntimeError:
            raise
        exposure = butler.get('calexp', dataId = {"visit": int(frameId), "ccd": ccdId})
        dCcdDims[ccdId] = exposure.getDimensions()
        dCcd[ccdId] = afwCG.cast_Ccd(exposure.getDetector())

    lZ = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.30]) + focusHeader - 0.3

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

    zemaxFilename = os.path.join(os.environ["OBS_SUBARU_DIR"], "hsc", "zemax_config%s_0.0.dat" % zemaxConfig)

    return testFocusDeterminationByFocusSweep(dSrc, dCcd, dCcdDims, lZ, zemaxFilename, frameId, config, plotFilename="tmp.png")

if __name__ == "__main__":
    frameId = int(sys.argv[1])
    rerun = "miyatake-miyatakesky-linear-lowdetectionthreshold-test"
    run(rerun, frameId)
