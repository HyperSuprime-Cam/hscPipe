# isrTasks for QA. 
# This module is based on ip_isr.isrTasks, modified for QA outputs

import lsst.afw.image       as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom  as cameraGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base import Struct
from lsst.ip.isr.ccdAssembler import CcdAssembler # equivelent to the above
from lsst.ip.isr import isrLib # equivelent to the above

## FH added for QA
import os, sys
import lsst.afw.math as afwMath
import lsst.ip.isr as ipIsr
import hsc.pipe.tasks.suprimecam as hscSuprimeCam
import hsc.pipe.tasks.isr as hscIsr
from hsc.pipe.tasks.qaSuprimeCamIsr import QaSuprimeCamIsr


#== FH changed for QA output
class QaFlatnessConfig(pexConfig.Config):
    meshX = pexConfig.Field(
        dtype = int,
        doc = 'Mesh size in X (pix) to calculate count statistics',
        default = 256,
        )
    meshY = pexConfig.Field(
        dtype = int,
        doc = 'Mesh size in Y (pix) to calculate count statistics',
        default = 256,
        )
    doClip = pexConfig.Field(
        dtype = bool,
        doc = 'Do we clip outliers in calculate count statistics?',
        default = True,
        )
    clipSigma = pexConfig.Field(
        dtype = float,
        doc = 'How many sigma is used to clip outliers in calculate count statistics?',
        default = 3.0,
        )
    nIter = pexConfig.Field(
        dtype = int,
        doc = 'How many times do we iterate clipping outliers in calculate count statistics?',
        default = 3,
        )

class QaDoWriteImageConfig(pexConfig.Config):
    doWriteOssImage = pexConfig.Field(
        dtype = bool,
        doc = 'Do we write overscan-subtracted image FITS?',
        default = True,
        )
    doWriteFltImage = pexConfig.Field(
        dtype = bool,
        doc = 'Do we write flatfielded image FITS?',
        default = True,
        )
    doDumpSnapshot = pexConfig.Field(
        dtype=bool,
        doc="Do we dump snapshot files?",
        default=True
        )

class QaConfig(pexConfig.Config):
    flatness = pexConfig.ConfigField(dtype=QaFlatnessConfig, doc="Qa.flatness")
    doWriteImage = pexConfig.ConfigField(dtype=QaDoWriteImageConfig, doc="Qa.DoWriteImage")
    camera = pexConfig.Field(dtype=str, doc="suprimecam or hsc", default='hsc')

#class QaIsrTaskConfig(ptProcessCcd.ProcessCcdConfig):
###class QaIsrTaskConfig(ipIsr.IsrTaskConfig):
###class QaIsrTaskConfig(hscIsr.HscIsrConfig):
class QaIsrTaskConfig(hscSuprimeCam.HamamatsuIsrTaskConfig):
    qa = pexConfig.ConfigField(dtype=QaConfig, doc="Qa configuration")

class QaSuprimeCamIsrTask(hscSuprimeCam.SuprimeCamIsrTask):

    ConfigClass = QaIsrTaskConfig

    def __init__(self, *args, **kwargs):
        ## We do not use methodList which is initialized in ip_isr.IsrTask.__init__() and
        ## other parent classes do not have specific __init__(), so I just call pipeBase.Task.__init__().
        ##super(QaSuprimeCamIsrTask, self).__init__(**kwargs)
        pipeBase.Task.__init__(self, *args, **kwargs)
        self.isr = QaSuprimeCamIsr()

    ## The run() function is copied from hsc.pipe.tasks.isr(=hscIsr).HscIsrTask.run() 
    ## with changing doWriteXxxImageQa() functions to accept <butler & dataId> in their arguments
    def run(self, exposure, calibSet, butler, dataId):
        exposure = self.doConversionForIsr(exposure, calibSet)
        if self.config.doSaturation:
            exposure = self.doSaturationDetection(exposure, calibSet)
        if self.config.doOverscan:
            exposure = self.doOverscanCorrectionQa(exposure, calibSet)
        if self.config.doVariance:
            # Ideally, this should be done after bias subtraction, but CCD assembly demands a variance plane
            exposure = self.doVariance(exposure, calibSet)
        ## FH: doWriteOssImageQa includes ccdAssembly which requires variance plance there, so
        ##     placed after doVariance
        if self.config.qa.doWriteImage.doWriteOssImage or self.config.qa.doWriteImage.doDumpSnapshot:
            self.doWriteOssImageQa(exposure, calibSet, butler, dataId)

        exposure = self.doCcdAssembly([exposure])

        if self.config.doBias:
            exposure = self.doBiasSubtraction(exposure, calibSet)
        if self.config.doDark:
            exposure = self.doDarkCorrection(exposure, calibSet)
        if self.config.doFlat:
            exposure = self.doFlatCorrectionQa(exposure, calibSet)
        if self.config.qa.doWriteImage.doWriteFltImage or self.config.qa.doWriteImage.doDumpSnapshot:
            self.doWriteFltImageQa(exposure, calibSet, butler, dataId)
        if self.config.qa.camera in ['suprimecam', 'suprime', 'sup', 'sc']:
            self.guider(exposure) ## to be compatible with hscSuprimeCam.SuprimeCamIsrTask.run()

        return Struct(postIsrExposure=exposure)

    def doCrosstalkCorrection(self, exposure, calibSet):
        pass

    ##== FH for QA output
    def doOverscanCorrectionQa(self, exposure, calibSet):
        fittype = self.config.overscanFitType
        polyorder = self.config.overscanPolyOrder

        metadata = exposure.getMetadata()
        osLevel = [None,None,None,None] # overscan values in each readout channel
        osSigma = [None,None,None,None]

        for amp in self._getAmplifiers(exposure):
            expImage = exposure.getMaskedImage().getImage()
            overscan = expImage.Factory(expImage, amp.getDiskBiasSec())
            exp = exposure.Factory(exposure, amp.getDiskDataSec())

            channelId = int(amp.getId().getSerial()) -1 # we should decide channelId is 0..3 or 1..4
            self.log.log(self.log.INFO, "QA overscan: channelId being processed: %d" % channelId)

            osLevel[channelId], osSigma[channelId] = self.isr.overscanCorrectionQa(exp.getMaskedImage(), overscan, fittype=fittype, polyorder=polyorder, imageFactory=afwImage.ImageF)

            self.log.log(self.log.INFO, "QA overscan: osLevel[%d]: %f" % (channelId, osLevel[channelId]))
            self.log.log(self.log.INFO, "QA overscan: osSigmal[%d]: %f" % (channelId, osSigma[channelId]))

        metadata.set('OSLEVEL1', osLevel[0])
        metadata.set('OSLEVEL2', osLevel[1])
        metadata.set('OSLEVEL3', osLevel[2])
        metadata.set('OSLEVEL4', osLevel[3])
        metadata.set('OSSIGMA1', osSigma[0])
        metadata.set('OSSIGMA2', osSigma[1])
        metadata.set('OSSIGMA3', osSigma[2])
        metadata.set('OSSIGMA4', osSigma[3])

        return exposure

##== FH for QA output
    def doWriteOssImageQa(self, exposure, calibSet, butler, dataId):
        trimmedExposure = self.doCcdAssembly([exposure])
        frameId = exposure.getMetadata().get('FRAMEID')	
        pathToSrcFile = butler.get('src_filename', dataId)[0]
        qaOutputDirName = os.path.dirname(pathToSrcFile)
        if os.path.exists(qaOutputDirName) is not True:
            os.makedirs(qaOutputDirName)
        else:
            pass
        if self.config.qa.doWriteImage.doWriteOssImage is True:
            outfile = qaOutputDirName + '/' + 'oss_'+frameId+'.fits'
            self.isr.writeFitsImageQa(trimmedExposure, outfile)
            self.log.log(self.log.INFO, "QA writing overscan-subtracted FITS image: %s" % outfile)

        #if True:
        if self.config.qa.doWriteImage.doDumpSnapshot is True:
            snapName = qaOutputDirName + '/' + 'oss_%s.png' % str(frameId)
            self.isr.writeSnapshotImageQa(trimmedExposure, snapName, format='png', width=500)
            self.log.log(self.log.INFO, "QA writing snapshot png of overscan-subtracted image: %s" % snapName)

        return exposure

    ##== FH for QA output
    def doWriteFltImageQa(self, exposure, calibSet, butler, dataId):
        frameId = exposure.getMetadata().get('FRAMEID')
        pathToSrcFile = butler.get('src_filename', dataId)[0]
        qaOutputDirName = os.path.dirname(pathToSrcFile)
        if os.path.exists(qaOutputDirName) is not True:
            os.makedirs(qaOutputDirName)
        else:
            pass
        if self.config.qa.doWriteImage.doWriteFltImage is True:
            outfile = qaOutputDirName + '/' + 'flt_'+frameId+'.fits'
            self.isr.writeFitsImageQa(exposure, outfile)
            self.log.log(self.log.INFO, "QA writing flatfielded FITS image: %s" % outfile)

        #if True:
        if self.config.qa.doWriteImage.doDumpSnapshot is True:
            snapName = qaOutputDirName + '/' + 'flt_%s.png' % str(frameId)
            self.isr.writeSnapshotImageQa(exposure, snapName, format='png', width=500)
            self.log.log(self.log.INFO, "QA writing snapshot png of flatfielded image: %s" % snapName)

        return exposure

    ##== FH for QA output
    def doFlatCorrectionQa(self, exposure, calibSet):
        flatfield = calibSet['flat']
        scalingtype = self.config.flatScalingType
        scalingvalue = self.config.flatScalingValue
        #### FH
        if False: # if self._display:
            import lsst.afw.display.ds9 as ds9

        for amp in self._getAmplifiers(exposure):
            exp, flat = self._getCalibration(exposure, flatfield, amp)
            self.isr.flatCorrection(exp.getMaskedImage(), flat.getMaskedImage(), scalingtype, scaling = scalingvalue)
            #### FH
            if False:
                ampId = amp.getId().getSerial()
                ds9.mtv(exp.getMaskedImage(), frame=2*ampId, title='exp-amp %1d' % ampId)
                ds9.mtv(flat.getMaskedImage(), frame=2*ampId+1, title='flat-amp %1d' % ampId)

        # Qa mesurement of skylevel of flatfielded image
        # - As background is not done here, I measure skyLevel and skySigma here.  
        self.doMeasureFlatnessImageQa(exposure)

        return exposure

    ##== FH for QA output
    def doMeasureFlatnessImageQa(self, exposure):

        metadata = exposure.getMetadata()
        ### exposure is already trimmed just before doFlatCorrection so no need to ccdAssemble here.
        if False:
            trimmedExposure = self.doCcdAssembly([exposure])
            trimmedImage = trimmedExposure.getMaskedImage().getImage()
        else:
            trimmedImage = exposure.getMaskedImage().getImage()

        clipSigma = 3.0; nIter = 3
        sctrl = afwMath.StatisticsControl(clipSigma, nIter)
        stats = afwMath.makeStatistics(trimmedImage, afwMath.MEDIAN|afwMath.STDEVCLIP, sctrl)
        skyLevel = stats.getValue(afwMath.MEDIAN)
        skySigma = stats.getValue(afwMath.STDEVCLIP)
        self.log.log(self.log.INFO, "QA skylevel and sigma: %f  %f" % (skyLevel, skySigma))
        metadata.set('SKYLEVEL', skyLevel)
        metadata.set('SKYSIGMA', skySigma)

        # Qa mesurement of flatness of flatfielded image
        print '***** self.config:'
        print self.config

        if True:
            configFlatness = self.config.qa.flatness
            meshX = configFlatness.meshX
            meshY = configFlatness.meshY
            doClip = configFlatness.doClip
            clipSigma = configFlatness.clipSigma
            nIter = configFlatness.nIter
        else:
            meshX = 256
            meshY = 256
            doClip = True
            clipSigma = 3.0
            nIter = 3

        (flatness, flatness_pp, flatness_min, flatness_max, flatness_rms, skyMedian, nX, nY) = \
                   self.isr.measureFlatnessImageQa(
            trimmedImage,
            meshX=meshX,
            meshY=meshY,
            doClip=doClip,
            clipSigma=clipSigma,
            nIter=nIter
            )

        self.log.log(self.log.INFO, "QA flatfield: measuring skylevels in %dx%d grids: %f" % (nX, nY, skyMedian))
        self.log.log(self.log.INFO, "QA flatfield: flatness in %dx%d grids - pp: %f rms: %f" % (nX, nY, flatness_pp, flatness_rms))

        metadata.set('FLATNESS_PP', flatness_pp)
        metadata.set('FLATNESS_RMS', flatness_rms)
        metadata.set('FLATNESS_NGRIDS', '%dx%d' % (nX, nY))
        metadata.set('FLATNESS_MESHX', meshX)
        metadata.set('FLATNESS_MESHY', meshY)
#    def doIlluminationCorrection(self, exposure, calibSet):
#        pass

