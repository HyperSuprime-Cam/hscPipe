from lsst.pex.config import Config, ConfigurableField, Field, ListField, ConfigField
from lsst.pipe.base import Task, Struct, ArgumentParser, TaskRunner
from lsst.pipe.base.timer import timeMethod
from hsc.pipe.base.pool import Pool
import lsst.afw.cameraGeom as cameraGeom
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.meas.algorithms as measAlg
import numpy
import lsst.obs.subaru.isr as hscIsr
import os
import math
import lsst.afw.display.ds9 as ds9


################################################################################
# LargeScaleBackgroundSubtractionTask
################################################################################

class LargeScaleBackgroundSubtractionConfig(Config):
    scale = Field(dtype=float, default=0.1, doc="scale for projection to focal plane")
    mask = ListField(doc="Mask planes to respect", dtype=str, default=["BAD", "NO_DATA"])
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    detectSigma = Field(dtype=float, default=5.0, doc="Detection PSF gaussian sigma")
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")


    def setDefaults(self):
        super(LargeScaleBackgroundSubtractionConfig, self).setDefaults()
        self.detection.reEstimateBackground = False
        self.detection.doFootprintBackground = True
        self.detection.footprintBackground.useApprox = True
        self.detection.footprintBackground.approxOrder = 3
        self.background.isNanSafe = True
        self.background.useApprox = False
        self.background.algorithm = 'AKIMA_SPLINE'
        self.background.binSize = 600


class LargeScaleBackgroundSubtractionTask(Task):
    _DefaultName = 'largeScaleBackgroundSubtraction'
    ConfigClass = LargeScaleBackgroundSubtractionConfig


    def __init__(self, *args, **kwargs):
        super(LargeScaleBackgroundSubtractionTask, self).__init__(*args, **kwargs)
        self.makeSubtask("detection")


    def run(self, pool, ccdRefList, butler):
        tiledImage = self._makeTiledImage(pool, ccdRefList, butler)
        background = self._makeLargeScaleBackground(tiledImage.maskedImage)
        pool.mapToPrevious(self._makeCcdBackground, ccdRefList, background, tiledImage.ccdLaytoutDict)


    def _makeLargeScaleBackground(self, maskedImage):
        self.log.info('making global background...')
        self.detection.detectFootprints(afwImage.ExposureF(maskedImage), sigma=self.config.detectSigma)
        background = measAlg.getBackground(maskedImage, self.config.background).getImageF()
        background.setXY0(maskedImage.getXY0())
        return background


    def _makeCcdBackground(self, cache, ccdRef, background, ccdLaytoutDict):
        ccd = ccdRef.dataId['ccd']
        self.log.info('subtracting global sky: ccd={}'.format(ccd))
        layout = ccdLaytoutDict[ccd]

        CX, CY = numpy.meshgrid( # ccd index
            numpy.arange(layout.width),
            numpy.arange(layout.height),
        )
        FX, FY = layout.focalPlaneCoordFromCcdCoord(CX, CY)

        # bilinear warp
        x0, y0 = background.getXY0()
        FXI, FXR = indexAndRatio(FX)
        FYI, FYR = indexAndRatio(FY)
        focalPlaneBackgroundArray = background.getArray()
        # I01 - I11
        #  |     |
        # I00 - I10
        I00 = imageIndexMap(focalPlaneBackgroundArray, FXI + x0,     FYI + y0)[0]
        I10 = imageIndexMap(focalPlaneBackgroundArray, FXI + x0 + 1, FYI + y0)[0]
        I01 = imageIndexMap(focalPlaneBackgroundArray, FXI + x0,     FYI + y0 + 1)[0]
        I11 = imageIndexMap(focalPlaneBackgroundArray, FXI + x0 + 1, FYI + y0 + 1)[0]
        ccdBackgroundArray = (
            (I00*(1 - FXR) + I10*FXR) * (1 - FYR) +
            (I01*(1 - FXR) + I11*FXR) * FYR
            )
        exposure = cache.postISRCCD
        array = exposure.getMaskedImage().getImage().getArray()
        array -= ccdBackgroundArray
        # exposure.writeFits('ccdbg/{}.fits'.format(ccd))


    def _makeTiledImage(self, pool, ccdRefList, butler):
        ccdLaytoutDict = CcdLayout.layoutDict(butler.mapper.camera, scale=self.config.scale)
        target = self._createTagetMaskedImage(ccdLaytoutDict)
        targetImageArray = target.getImage().getArray()
        ok = numpy.zeros_like(targetImageArray, dtype=bool)
        x0, y0 = target.getXY0()

        warpedList = pool.mapToPrevious(self._warpToFocalPlane, ccdRefList, ccdLaytoutDict)

        for i, warpResult in enumerate(warpedList):
            if ccdRefList[i].dataId['ccd'] == 9:
                continue
            yRagne = warpResult.focalYIndexRange
            xRagne = warpResult.focalXIndexRange
            targetImageArray[
                yRagne[0] - y0 : yRagne[1] - y0,
                xRagne[0] - x0 : xRagne[1] - x0 ] += warpResult.warpedImageArray
            ok[
                yRagne[0] - y0 : yRagne[1] - y0,
                xRagne[0] - x0 : xRagne[1] - x0 ] |=  ~ warpResult.noData

        noDataBit = 1 << target.getMask().addMaskPlane('NO_DATA')
        target.getMask().getArray()[~ ok] = noDataBit

        return Struct(
            maskedImage=target,
            ccdLaytoutDict=ccdLaytoutDict,
        )


    def _warpToFocalPlane(self, cache, ccdRef, ccdLaytoutDict):
        self.log.info('warping ccd={}...'.format(ccdRef.dataId['ccd']))
        layout = ccdLaytoutDict[ccdRef.dataId['ccd']]
        bbox = bboxDouble2Int(layout.bboxOnFocalPlane())
        focalXIndexRange = (bbox.getMinX(), bbox.getMaxX() + 1)
        focalYIndexRange = (bbox.getMinY(), bbox.getMaxY() + 1)

        FX, FY = numpy.meshgrid( # focal plane coord
            numpy.arange(*focalXIndexRange),
            numpy.arange(*focalYIndexRange),
            )
        CX, CY = layout.ccdCoordFromFocalPlaneCoord(FX, FY) # ccd coord index
        CX = npFloorInt(CX) # ccd index
        CY = npFloorInt(CY)

        srcExposure = cache.postISRCCD

        srcImageArray = srcExposure.getMaskedImage().getImage().getArray()
        warpedImageArray, noData = imageIndexMap(srcImageArray, CX, CY)
        warpedImageArray[noData] = 0.0

        srcMaskArray = srcExposure.getMaskedImage().getMask().getArray()
        warpedMaskArray, noData = imageIndexMap(srcMaskArray, CX, CY)
        noData = numpy.logical_or(noData, warpedMaskArray & self._errorBit(srcExposure))

        return Struct(
            warpedImageArray=warpedImageArray,
            noData=noData,
            focalXIndexRange=focalXIndexRange,
            focalYIndexRange=focalYIndexRange,
        )


    def _createTagetMaskedImage(self, ccdLaytoutDict):
        corners = []
        for layout in ccdLaytoutDict.values():
            bbox = layout.bboxOnFocalPlane()
            corners.append([bbox.getMinX(), bbox.getMinY()])
            corners.append([bbox.getMaxX(), bbox.getMaxY()])
        corners = numpy.array(corners)
        minX, minY = corners.min(axis=0)
        maxX, maxY = corners.max(axis=0)
        bbox = afwGeom.Box2D(
            afwGeom.Point2D(minX, minY),
            afwGeom.Point2D(maxX, maxY))

        target = afwImage.MaskedImageF(bboxDouble2Int(bbox))
        return target


    def _errorBit(self, exposure):
        maskDict = dict(exposure.getMaskedImage().getMask().getMaskPlaneDict())
        return sum((1 << maskDict[plane] for plane in self.config.mask), 0)



################################################################################
# SkyPatternSubtractionTask
################################################################################

class SkyPatternSubtractionConfig(Config):
    mask = ListField(doc="Mask planes to respect", dtype=str, default=["BAD", "NO_DATA", "DETECTED"])
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    detectSigma = Field(dtype=float, default=5.0, doc="Detection PSF gaussian sigma")

    def setDefaults(self):
        super(SkyPatternSubtractionConfig, self).setDefaults()
        self.detection.reEstimateBackground = False
        self.detection.doFootprintBackground = True
        self.detection.footprintBackground.useApprox = True
        self.detection.footprintBackground.approxOrder = 3


class SkyPatternSubtractionTask(Task):
    _DefaultName = 'skyPatternSubtraction'
    ConfigClass = SkyPatternSubtractionConfig

    def __init__(self, *args, **kwargs):
        super(SkyPatternSubtractionTask, self).__init__(*args, **kwargs)
        self.makeSubtask('detection')


    def run(self, pool, ccdRefList, butler, skyPatternDir):
        median = self._computeMedian(pool, ccdRefList)
        pool.mapToPrevious(self._subtractSkyPattern, ccdRefList, median, skyPatternDir)


    def _computeMedian(self, pool, ccdRefList):
        self.log.info('getting median')
        arrays = pool.mapToPrevious(self._getSafeValues, ccdRefList)
        return numpy.median(numpy.concatenate(arrays))


    def _getSafeValues(self, cache, ccdRef):
        ccd = ccdRef.dataId['ccd']
        self.log.info('getting safe values: ccd = {}'.format(ccd))
        if ccd == 9:
            return numpy.array([])
        postISRCCD = cache.postISRCCD.__class__(cache.postISRCCD, True) # copy
        self.detection.detectFootprints(postISRCCD, sigma=self.config.detectSigma)
        mi = postISRCCD.getMaskedImage()
        imageArray = mi.getImage().getArray()
        maskArray = mi.getMask().getArray()
        return imageArray[maskArray & self._errorBit(postISRCCD) == 0]


    def _subtractSkyPattern(self, cache, ccdRef, median, skyPatternDir):
        ccd = ccdRef.dataId['ccd']
        self.log.info('subtracting sky pattern: ccd = {}'.format(ccd))
        skyPattern = afwImage.ExposureF('{}/{}.fits'.format(skyPatternDir, ccd))
        skyPatternMi = skyPattern.getMaskedImage()
        self.log.info('median={}'.format(median))
        skyPatternMi *= float(median)
        postISRCCDMi = cache.postISRCCD.getMaskedImage()
        postISRCCDMi -= skyPatternMi


    def _errorBit(self, exposure):
        maskDict = dict(exposure.getMaskedImage().getMask().getMaskPlaneDict())
        return sum((1 << maskDict[plane] for plane in self.config.mask), 0)


################################################################################
# HscIsrWithGlobalSkySubtractionTask
################################################################################

class HscIsrWithGlobalSkySubtractionConfig(Config):
    isr = ConfigurableField(target=hscIsr.SubaruIsrTask, doc="ISR configuration")
    skyPatternSubtraction = ConfigurableField(target=SkyPatternSubtractionTask, doc='sky pattern sub')
    largeScaleBackgroundSubtraction = ConfigurableField(target=LargeScaleBackgroundSubtractionTask, doc='large scale bg sub')
    doSkyPatternSubtraction = Field(dtype=bool, default=True, doc="do sky pattern sub?")
    doLargeScaleBackgroundSubtraction = Field(dtype=bool, default=True, doc="do large scale background sub?")
    doWrite = Field(dtype=bool, default=True, doc="do write?")

    def setDefaults(self):
        super(HscIsrWithGlobalSkySubtractionConfig, self).setDefaults()
        self.load(os.environ['OBS_SUBARU_DIR'] + '/config/hsc/isr.py')


class HscIsrWithGlobalSkySubtractionTask(Task):
    _DefaultName = 'hscIsrWithGlobalSkySubtraction'
    ConfigClass = HscIsrWithGlobalSkySubtractionConfig


    def __init__(self, *args, **kwargs):
        super(HscIsrWithGlobalSkySubtractionTask, self).__init__(*args, **kwargs)
        self.makeSubtask("isr")
        self.makeSubtask("largeScaleBackgroundSubtraction")
        self.makeSubtask("skyPatternSubtraction")


    def run(self, pool, ccdRefList, butler, skyPatternDir):
        if self.config.doSkyPatternSubtraction:
            assert skyPatternDir
        pool.map(self._isr, ccdRefList)
        # pool.mapToPrevious(self._debugWrite, ccdRefList, 'before-skypatternSub')
        if self.config.doSkyPatternSubtraction:
            self.skyPatternSubtraction.run(pool, ccdRefList, butler, skyPatternDir)
            # pool.mapToPrevious(self._debugWrite, ccdRefList, 'skypatternSub')
        if self.config.doLargeScaleBackgroundSubtraction:
            self.largeScaleBackgroundSubtraction.run(pool, ccdRefList, butler)
            # pool.mapToPrevious(self._debugWrite, ccdRefList, 'largeScaleBGSub')
        if self.config.doWrite:
            pool.mapToPrevious(self._writePostISRCCD, ccdRefList)
        pool.cacheClear()


    def _isr(self, cache, ccdRef):
        # if ccdRef.datasetExists('postISRCCD'):
        #     cache.postISRCCD = ccdRef.get('postISRCCD')
        # else:
        #     cache.postISRCCD = self.isr.run(ccdRef).exposure
        #     ccdRef.put(cache.postISRCCD, 'postISRCCD')
        cache.postISRCCD = self.isr.run(ccdRef).exposure
        ccdRef.put(cache.postISRCCD, 'postISRCCD')


    def _writePostISRCCD(self, cache, ccdRef):
        ccdRef.put(cache.postISRCCD, 'postISRCCD')


    # def _debugWrite(self, cache, ccdRef, outDir):
    #     try:
    #         os.makedirs(outDir)
    #     except:
    #         pass
    #     cache.postISRCCD.writeFits('{}/{}.fits'.format(outDir, ccdRef.dataId['ccd']))


################################################################################
# utilities
################################################################################

def bboxDouble2Int(bbox):
    return afwGeom.Box2I(
        afwGeom.Point2I(floorInt(bbox.getMinX()),     floorInt(bbox.getMinY())),
        afwGeom.Point2I(floorInt(bbox.getMaxX() + 1), floorInt(bbox.getMaxY() + 1)),
        )


def floorInt(x):
    return int(math.floor(x))


def npFloorInt(x):
    return numpy.array(numpy.floor(x), dtype=int)


def indexAndRatio(x):
    floor = numpy.floor(x)
    r = x - floor
    return (numpy.array(floor, dtype=int), r)


def imageIndexMap(image, X, Y):
    h, w = image.shape
    noData = numpy.logical_or(
        numpy.logical_or(Y < 0, Y >= image.shape[0]),
        numpy.logical_or(X < 0, X >= image.shape[1]),
    )
    I = (Y * w + X).flatten()
    I %= image.size
    warped = image.flatten()[I].reshape(X.shape)
    return (warped, noData)


class CcdLayout(object):
    def __init__(self, ccdGeom, scale=1):
        self.id = ccdGeom.getId().getSerial()
        self.width, self.height = ccdGeom.getAllPixels(True).getDimensions()
        o = ccdGeom.getPositionFromPixel(afwGeom.Point2D(0, 0)).getPixels(1.) # focal plane position
        x = ccdGeom.getPositionFromPixel(afwGeom.Point2D(1, 0)).getPixels(1.) - o
        y = ccdGeom.getPositionFromPixel(afwGeom.Point2D(0, 1)).getPixels(1.) - o
        self.affine = AffineTransform(scale * numpy.array([[x[0], y[0]], [x[1], y[1]]]), scale * numpy.array(o))


    def focalPlaneCoordFromCcdCoord(self, x, y):
        return self.affine.apply(x, y)


    def ccdCoordFromFocalPlaneCoord(self, x, y):
        return self.affine.inverse().apply(x, y)


    def bboxOnFocalPlane(self):
        corners = numpy.array([
            self.focalPlaneCoordFromCcdCoord(x, y)
            for x, y in ((0, 0), (self.width, 0), (self.width, self.height), (0, self.height))
        ])
        minX, minY = corners.min(axis=0)
        maxX, maxY = corners.max(axis=0)
        return afwGeom.Box2D(
            afwGeom.Point2D(minX, minY),
            afwGeom.Point2D(maxX, maxY))


    @classmethod
    def layoutDict(cls, camera, scale=1):
        layoutDict = {}
        for raft in camera:
            raft = cameraGeom.cast_Raft(raft)
            for ccdGeom in raft:
                ccdGeom = cameraGeom.cast_Ccd(ccdGeom)
                ccdGeom.setTrimmed(True)
                ccdLayout = CcdLayout(ccdGeom, scale)
                layoutDict[ccdGeom.getId().getSerial()] = ccdLayout
        return layoutDict


class AffineTransform(object):
    def __init__(self, a, b):
        self.a = a
        self.b = b


    def apply(self, x, y):
        # A x + b
        a = self.a
        b = self.b
        return (
            a[0][0]*x + a[0][1]*y + b[0],
            a[1][0]*x + a[1][1]*y + b[1],)


    def inverse(self):
        c = numpy.linalg.inv(self.a)
        return AffineTransform(c, - c.dot(self.b))
