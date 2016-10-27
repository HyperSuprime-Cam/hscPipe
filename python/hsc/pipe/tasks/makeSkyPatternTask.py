from lsst.pex.config import Config, ConfigField, ConfigurableField, Field, ListField
from lsst.pipe.base import Task, Struct, ArgumentParser, TaskRunner
import lsst.afw.detection as afwDet
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.obs.subaru.isr as hscIsr
import numpy
import os
from hsc.pipe.base.parallel import BatchPoolTask
from hsc.pipe.base.pool import Pool, NODE, Debugger as PoolDebugger
import hsc.pipe.base.butler as hscButler


class SmoothConfig(Config):
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background configuration")

    def setDefaults(self):
        super(SmoothConfig, self).setDefaults()
        self.detection.thresholdValue=3.0
        self.detection.reEstimateBackground = False
        self.detection.doFootprintBackground = True
        self.detection.footprintBackground.useApprox = True
        self.background.useApprox = False
        self.background.binSize = 256


class SmoothTask(Task):
    _DefaultName = 'smooth'
    ConfigClass = SmoothConfig

    def __init__(self, *args, **kwargs):
        super(SmoothTask, self).__init__(*args, **kwargs)
        self.makeSubtask("detection")

    def run(self, exp):
        tmpExp = exp.__class__(exp, True) # copy
        self.detection.detectFootprints(tmpExp, sigma=5.)
        image = exp.getMaskedImage().getImage()
        image <<= measAlg.getBackground(tmpExp.getMaskedImage(), self.config.background).getImageF()


class MakeSkyPatternConfig(Config):
    isr = ConfigurableField(target=hscIsr.SubaruIsrTask, doc="ISR configuration")
    doDetection = Field(doc="do detection?", dtype=bool, default=True)
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Detection configuration")
    detectSigma = Field(dtype=float, default=5.0, doc="Detection PSF gaussian sigma")
    mask = ListField(doc="Mask planes to respect", dtype=str, default=["BAD", "SAT", "DETECTED", "INTRP", "NO_DATA"])
    combine = Field(doc="Statistic to use for combination (from lsst.afw.math)", dtype=int, default=afwMath.MEDIAN)
    clip = Field(doc="Clipping threshold for combination", dtype=float, default=3.0)
    iter = Field(doc="Clipping iterations for combination", dtype=int, default=3)
    rows = Field(doc="Number of rows to read at a time", dtype=int, default=512)
    doSmooth = Field(doc="do smooth?", dtype=bool, default=True)
    smooth = ConfigurableField(target=SmoothTask, doc="smooth configuration")

    def setDefaults(self):
        super(MakeSkyPatternConfig, self).setDefaults()
        self.load(os.environ['OBS_SUBARU_DIR'] + '/config/hsc/isr.py')
        self.detection.reEstimateBackground = False
        self.detection.doFootprintBackground = True
        self.detection.footprintBackground.useApprox = True


class MakeSkyPatternTaskRunner(TaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [dict(refList=parsedCmd.id.refList, parsedCmd=parsedCmd)]

    def __call__(self, args):
        task = self.TaskClass(config=self.config, log=self.log)
        task.run(**args)


class MakeSkyPatternTask(BatchPoolTask):
    _DefaultName = "makeSkyPattern"
    RunnerClass = MakeSkyPatternTaskRunner
    ConfigClass = MakeSkyPatternConfig


    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", help="data ID, e.g. --id visit=12345")
        parser.add_argument('--sky-pattern-dir', '-O', required=True)
        return parser


    def __init__(self, *args, **kwargs):
        super(MakeSkyPatternTask, self).__init__(*args, **kwargs)
        self.makeSubtask("detection")
        self.makeSubtask("isr")
        self.makeSubtask("smooth")


    def run(self, refList, parsedCmd):
        self.visitList = sorted(set(ref.dataId['visit'] for ref in refList))
        self.ccdList = sorted(set(ref.dataId['ccd'] for ref in refList))
        assert len(self.visitList) * len(self.ccdList) == len(refList)

        self.pool = Pool()
        self.pool.storeSet(
            butler=parsedCmd.butler,
            outputDir=parsedCmd.sky_pattern_dir,
        )

        # isr & mask objects
        self.log.info('doing isr...')
        self._isr(refList)

        # determine scales
        self.log.info('determining scales...')
        visitMedianDict = self._getVisitMedian(refList)
        self.log.info('visitMedianDict: {}'.format(visitMedianDict))

        # combine visits
        self._combineAndWrite(refList, visitMedianDict)


    def _combineAndWrite(self, refList, visitMedianDict):
        # batch by ccd
        ccdData = {}
        for ccd in self.ccdList:
            ccdData[ccd] = [ref for ref in refList if ref.dataId['ccd'] == ccd ]
        self.pool.map(self._combineAndWriteSingle, ccdData.values(), visitMedianDict)


    def _combineAndWriteSingle(self, cache, refList, visitMedianDict):
        ccd = refList[0].dataId['ccd']
        self.log.info('stacking ccd={}...'.format(ccd))

        maskVal = 0
        for mask in self.config.mask:
            maskVal |= afwImage.MaskU.getPlaneBitMask(mask)
        stats = afwMath.StatisticsControl(self.config.clip, self.config.iter, maskVal)

        templateExposure = refList[0].get('postISRCCD')
        templateExposure.getMaskedImage().getVariance().set(0.)
        width, height = templateExposure.getDimensions()
        combined = afwImage.ImageF(width, height)
        numImages = len(refList)
        imageList = [None] * numImages

        for start in range(0, height, self.config.rows):
            rows = min(self.config.rows, height - start)
            box = afwGeom.Box2I(afwGeom.Point2I(0, start), afwGeom.Extent2I(width, rows))
            subCombined = combined.Factory(combined, box)
            for i, ref in enumerate(refList):
                exposure = ref.get("postISRCCD_sub", bbox=box)
                mi = exposure.getMaskedImage()
                mi /= visitMedianDict[ref.dataId['visit']]
                imageList[i] = mi
            subCombined <<= afwMath.statisticsStack(
                afwImage.vectorMaskedImageF(imageList),
                self.config.combine,
                stats).getImage()

        self._makeCombinedExposure(templateExposure, combined)
        self._writeCombined(templateExposure, ccd, cache.outputDir)


    def _makeCombinedExposure(self, template, combined):
        image = template.getMaskedImage().getImage()
        image <<= combined
        imageArray = image.getArray()
        mask = template.getMaskedImage().getMask()
        maskArray = mask.getArray()

        if self.config.doSmooth:
            self.smooth.run(template)
        else:
            bad = numpy.isnan(imageArray)
            maskArray.fill(0)
            maskArray[bad] = 1 << mask.addMaskPlane('BAD')
            imageArray[bad] = numpy.median(imageArray[numpy.logical_not(bad)])


    def _writeCombined(self, combined, ccd, outputDir):
        try:
            os.makedirs(outputDir)
        except:
            pass
        combined.writeFits('{}/{}.fits'.format(outputDir, ccd))


    def _isr(self, refList):
        self.pool.map(self._isrSingle, refList)


    def _isrSingle(self, cache, ref):
        # isr
        exposure = self.isr.run(ref).exposure
        # mask objects
        if self.config.doDetection:
            self.detection.detectFootprints(exposure, sigma=self.config.detectSigma)
        # write to disk
        ref.put(exposure, 'postISRCCD')


    def _getVisitMedian(self, refList):
        # batch by visit
        visitData = {}
        for visit in self.visitList:
            visitData[visit] = [ref for ref in refList if ref.dataId['visit'] == visit]
        medianList = self.pool.map(self._getVisitMedianSingle, visitData.values())
        return dict(zip(self.visitList, medianList))


    def _getVisitMedianSingle(self, cache, refList):
        arrays = []
        self.log.info('determining median of visit=%d...' % refList[0].dataId['visit'])
        for ref in refList:
            if ref.dataId['ccd'] == 9:
                continue
            exp = ref.get('postISRCCD')
            imageArray = exp.getMaskedImage().getImage().getArray()
            maskArray = exp.getMaskedImage().getMask().getArray()
            arrays.append(imageArray[numpy.logical_or(
                numpy.isnan(imageArray),
                maskArray & self._errorBit(exp)
            )])
        return float(numpy.median(numpy.concatenate(arrays)))


    def _errorBit(self, exposure):
        maskDict = dict(exposure.getMaskedImage().getMask().getMaskPlaneDict())
        return sum((1 << maskDict[plane] for plane in self.config.mask), 0)


    def _getConfigName(self):
        return None


if __name__ == '__main__':
    MakeSkyPatternTask.parseAndRun()
