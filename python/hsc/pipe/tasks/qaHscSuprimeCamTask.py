# isrTasks for QA. 
# This module is based on ip_isr.isrTasks, modified for QA outputs

import lsst.afw.image       as afwImage
import lsst.meas.algorithms as measAlg
import lsst.afw.cameraGeom  as cameraGeom
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
#from .isr import Linearization
#from .isr import Isr                 # commented  for QA inherit
#####from lsst.ip.isr.isr import Isr  # equivelent to the above
#from .ccdAssembler import CcdAssembler  
from lsst.ip.isr.ccdAssembler import CcdAssembler # equivelent to the above
#from . import isrLib
from lsst.ip.isr import isrLib # equivelent to the above

## FH added for QA
import lsst.afw.math as afwMath
import lsst.ip.isr as ipIsr
import hsc.pipe.tasks.suprimecam as hscSuprimeCam
from hsc.pipe.tasks.qaHscSuprimeCamIsr import qaSuprimeCamIsr


class IsrTaskConfig(pexConfig.Config):
    doWrite = pexConfig.Field(dtype=bool, doc="Write output?", default=True)
    fwhm = pexConfig.Field(
        dtype = float,
        doc = "FWHM of PSF (arcsec)",
        default = 1.0,
    )
    #This is needed for both the detection and correction aspects
    saturatedMaskName = pexConfig.Field(
        dtype = str,
        doc = "Name of mask plane to use in saturation detection",
        default = "SAT",
    )
    flatScalingType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for scaling the flat on the fly.",
        default = 'USER',
        allowed = {"USER": "User defined scaling",
            "MEAN": "Scale by the inverse of the mean",
            "MEDIAN": "Scale by the inverse of the median",
        },
    )
    flatScalingValue = pexConfig.Field(
        dtype = float,
        doc = "If scaling type is USER, a value for the scaling must be provided",
        default = 1.0,
    )
    overscanFitType = pexConfig.ChoiceField(
        dtype = str,
        doc = "The method for fitting the overscan bias level.",
        default = 'MEDIAN',
        allowed = {"POLY": "Fit polynomial to the longest axis of the overscan region",
            "MEAN": "Correct using the mean of the overscan region",
            "MEDIAN": "Correct using the median of the overscan region",
        },
    )
    overscanPolyOrder = pexConfig.Field(
        dtype = int,
        doc = "Order of polynomial to fit if overscan fit type is POLY",
        default = 1,
    )
    growSaturationFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the saturation footprints",
        default = 1,
    )
    growDefectFootprintSize = pexConfig.Field(
        dtype = int,
        doc = "Number of pixels by which to grow the defect (bad and nan) footprints",
        default = 1,
    )
    setGainAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "update exposure metadata in the assembled ccd to reflect the effective gain of the assembled chip",
        default = True,
    )
    keysToRemoveFromAssembledCcd = pexConfig.ListField(
        dtype = str,
        doc = "fields to remove from the metadata of the assembled ccd.",
        default = [],
    )
    reNormAssembledCcd = pexConfig.Field(
        dtype = bool,
        doc = "renormalize the assembled chips to have unity gain.  False if setGain is False",
        default = True,
    )
    methodList = pexConfig.ListField(
        dtype = str,   
        doc = "The list of ISR corrections to apply in the order they should be applied",
        default = ["doConversionForIsr", "doSaturationDetection", "doOverscanCorrection", "doBiasSubtraction", "doVariance", "doDarkCorrection", "doFlatCorrection"],
    )
    
#class IsrTask(pipeBase.Task):
class qaSuprimeCamIsrTask(hscSuprimeCam.SuprimeCamIsrTask):
    ConfigClass = IsrTaskConfig
    def __init__(self, *args, **kwargs):
        pipeBase.Task.__init__(self, *args, **kwargs)
        ## FH changed for QA output
        ##self.isr = Isr()
        self.isr = qaSuprimeCamIsr()
        self.methodList = []
        for methodname in self.config.methodList:
            self.methodList.append(getattr(self, methodname))

    def run(self, exposure, calibSet):
        """Do instrument signature removal on an exposure: saturation, bias, overscan, dark, flat, fringe correction

        @param exposure Apply ISR to this Exposure
        @param calibSet Dictionary of calibration products (bias/zero, dark, flat, fringe, linearization information)
        @return a pipeBase.Struct with fields:
        - postIsrExposure: the exposure after application of ISR
        """

        #The ISR routines operate in place.  A copy of the original exposure
        #will be made and the reduced exposure will be returned.
        workingExposure = exposure.Factory(exposure, True)
        for m in self.methodList:
            workingExposure = m(workingExposure, calibSet)
        return pipeBase.Struct(postIsrExposure=workingExposure)

    def runButler(self, butler, dataid):
        """Run the ISR given a butler
        @param butler Butler describing the data repository
        @param dataid A data identifier of the amp to process
        @return a pieBase.Struct see self.run for returned fields
        """
        calibSet = self.makeCalibDict(butler, dataid)
        output = self.run(butler.get("raw", dataid), calibSet)
        if self.config.doWrite:
            butler.put(output.postIsrExposure, "postISR", dataId=dataid)
        return output
        
    def makeCalibDict(self, butler, dataId):
        ret = {}
        required = {"doBiasSubtraction": "bias",
                    "doDarkCorrection": "dark",
                    "doFlatCorrection": "flat",
                    }
        for method in required.keys():
            if method in self.config.methodList:
                calib = required[method]
                ret[calib] = butler.get(calib, dataId)
        return ret

    def doConversionForIsr(self, exposure, calibSet):
        if not isinstance(exposure, afwImage.ExposureU):
            raise Exception("ipIsr.convertImageForIsr: Expecting Uint16 image. Got\
                %s."%(exposure.__repr__()))

        newexposure = exposure.convertF()
        mi = newexposure.getMaskedImage()
        var = afwImage.ImageF(mi.getBBox(afwImage.PARENT))
        mask = afwImage.MaskU(mi.getBBox(afwImage.PARENT))
        mask.set(0)
        newexposure.setMaskedImage(afwImage.MaskedImageF(mi.getImage(), mask, var))
        return newexposure

    def doVariance(self, exposure, calibSet):
        for amp in self._getAmplifiers(exposure):
            exp = exposure.Factory(exposure, amp.getDiskDataSec())
            self.isr.updateVariance(exp.getMaskedImage(), amp.getElectronicParams().getGain())
        return exposure

    def doCrosstalkCorrection(self, exposure, calibSet):
        pass

    def doSaturationCorrection(self, exposure, calibSet):
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        maskname = self.config.saturatedMaskName
        for amp in self._getAmplifiers(exposure):
            ep = amp.getElectronicParams()
            satvalue = ep.getSaturationLevel()
            exp = exposure.Factory(exposure, amp.getDiskDataSec())
            self.isr.saturationCorrection(exp.getMaskedImage(), satvalue, fwhm, growFootprints=grow, maskName=maskname)
        return exposure

    def doSaturationDetection(self, exposure, calibSet):
        for amp in self._getAmplifiers(exposure):
            datasec = amp.getDiskDataSec()
            exp = exposure.Factory(exposure, datasec)
            ep = amp.getElectronicParams()
            satvalue = ep.getSaturationLevel()
            maskname = self.config.saturatedMaskName
            self.isr.makeThresholdMask(exp.getMaskedImage(), satvalue, growFootprints=0, maskName=maskname)
        return exposure

    def doSaturationInterpolation(self, exposure, calibSet):
        #Don't loop over amps since saturation can cross amp boundaries
        maskname = self.config.saturatedMaskName
        fwhm = self.config.fwhm
        grow = self.config.growSaturationFootprintSize
        self.isr.interpolateFromMask(exposure.getMaskedImage(), fwhm, growFootprints=grow, maskName=maskname)
        return exposure
    
    def doMaskAndInterpDefect(self, exposure, calibSet):
        #Don't loop over amps since defects could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        defectBaseList = cameraGeom.cast_Ccd(exposure.getDetector()).getDefects()
        defectList = measAlg.DefectListT()
        #mask bad pixels in the camera class
        #create master list of defects and add those from the camera class
        for d in defectBaseList:
            bbox = d.getBBox()
            nd = measAlg.Defect(bbox)
            defectList.append(nd)
        self.isr.maskPixelsFromDefectList(exposure.getMaskedImage(), defectList, maskName='BAD')
        defectList = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='BAD', growFootprints=grow)
        self.isr.interpolateDefectList(exposure.getMaskedImage(), defectList, fwhm)
        return exposure

    def doMaskAndInterpNan(self, exposure, calibSet):
        #Don't loop over amps since nans could cross amp boundaries
        fwhm = self.config.fwhm
        grow = self.config.growDefectFootprintSize
        #find unmasked bad pixels and mask them
        exposure.getMaskedImage().getMask().addMaskPlane("UNMASKEDNAN") 
        unc = isrLib.UnmaskedNanCounterF()
        unc.apply(exposure.getMaskedImage())
        nnans = unc.getNpix()
        expmeta = exposure.getMetadata()
        expmeta.set("NUMNANS", nnans)
        if not nnans == 0:
		raise RuntimeError("There were %i unmasked NaNs"%(nnans))
        #get footprints of bad pixels not in the camera class
        undefects = self.isr.getDefectListFromMask(exposure.getMaskedImage(), maskName='UNMASKEDNAN', growFootprints=grow)
        #interpolate all bad pixels
        self.isr.interpolateDefectList(exposure.getMaskedImage(), undefects, fwhm)
        return exposure

    '''
    def doLinearization(self, exposure, calibSet):
        linearizer = Linearization(calibSet['linearityFile'])
        linearizer.apply(exposure)
        return exposure
    '''
    
#    def doOverscanCorrection(self, exposure, calibSet):
#        fittype = self.config.overscanFitType
#        polyorder = self.config.overscanPolyOrder
#        for amp in self._getAmplifiers(exposure):
#            expImage = exposure.getMaskedImage().getImage()
#            overscan = expImage.Factory(expImage, amp.getDiskBiasSec())
#            exp = exposure.Factory(exposure, amp.getDiskDataSec())
#            self.isr.overscanCorrection(exp.getMaskedImage(), overscan, fittype=fittype, polyorder=polyorder,
#                                        imageFactory=afwImage.ImageF)
#        return exposure

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

    def doBiasSubtraction(self, exposure, calibSet):
        biasExposure = calibSet['bias']
        for amp in self._getAmplifiers(exposure):
            exp, bias = self._getCalibration(exposure, biasExposure, amp)
            self.isr.biasCorrection(exp.getMaskedImage(), bias.getMaskedImage())
        
        return exposure 

    def doDarkCorrection(self, exposure, calibSet):
        darkexposure = calibSet['dark']
        darkscaling = darkexposure.getCalib().getExptime()
        expscaling = exposure.getCalib().getExptime()
        
        for amp in self._getAmplifiers(exposure):
            exp, dark = self._getCalibration(exposure, darkexposure, amp)
            self.isr.darkCorrection(exp.getMaskedImage(), dark.getMaskedImage(), expscaling, darkscaling)
        return exposure

    def doFringeCorrection(self, exposure, calibSet):
        pass

    def _getAmplifiers(self, exposure):
        """Return list of all amplifiers in an Exposure"""
        amp = cameraGeom.cast_Amp(exposure.getDetector())
        if amp is not None:
            return [amp]
        ccd = cameraGeom.cast_Ccd(exposure.getDetector())
        assert ccd is not None
        return [cameraGeom.cast_Amp(a) for a in ccd]

    def _getCalibration(self, exposure, calibration, amp):
        """Get a suitably-sized calibration exposure"""
        exp = exposure
        calib = calibration
        if exp.getDimensions() != calib.getDimensions():
            # Try just the exposure's pixels of interest
            try:
                exp = exp.Factory(exp, amp.getDiskDataSec()) # Exposure not trimmed or assembled
            except:
                pass
        if exp.getDimensions() != calib.getDimensions():
            # Try just the calibration's pixels of interest
            try:
                calib = calib.Factory(calib, amp.getDataSec(True)) # Calib is likely trimmed and assembled
            except:
                pass
        if exp.getDimensions() != calib.getDimensions():
            raise RuntimeError("Dimensions for exposure (%s) and calibration (%s) don't match" % \
                               (exposure.getDimensions(), calibration.getDimensions()))
        return exp, calib

    def doFlatCorrection(self, exposure, calibSet):
        flatfield = calibSet['flat']
        scalingtype = self.config.flatScalingType
        scalingvalue = self.config.flatScalingValue

        for amp in self._getAmplifiers(exposure):
            exp, flat = self._getCalibration(exposure, flatfield, amp)
            self.isr.flatCorrection(exp.getMaskedImage(), flat.getMaskedImage(), scalingtype, scaling = scalingvalue)   
        return exposure

    def doIlluminationCorrection(self, exposure, calibSet):
        pass

    def doCcdAssembly(self, exposureList):
        renorm = self.config.reNormAssembledCcd
        setgain = self.config.setGainAssembledCcd
        k2rm = self.config.keysToRemoveFromAssembledCcd
        assembler = CcdAssembler(exposureList, reNorm=renorm, setGain=setgain, keysToRemove=k2rm)
        return assembler.assembleCcd()
