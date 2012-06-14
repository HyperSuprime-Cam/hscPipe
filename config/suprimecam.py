root.doWriteIsr = True
##root.doWriteIsr = False
root.isr.methodList=["doConversionForIsr", "doSaturationDetection",
#                     "doOverscanCorrection", "doVariance", "doFlatCorrection"]
                     "doOverscanCorrectionQa", "doVariance", "doWriteOssImageQa", "doFlatCorrectionQa", "doWriteFltImageQa"]
##root.isr.doWrite = False
##root.isr.doWrite = True

root.calibrate.repair.doCosmicRay = True
root.calibrate.repair.cosmicray.nCrPixelMax = 1000000
root.calibrate.background.binSize = 256

# PSF determination
root.calibrate.measurePsf.starSelector.name = "secondMoment"
root.calibrate.measurePsf.psfDeterminer.name = "pca"
root.calibrate.measurePsf.starSelector["secondMoment"].clumpNSigma = 2.0
root.calibrate.measurePsf.psfDeterminer["pca"].nEigenComponents = 4
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSize = 10
root.calibrate.measurePsf.psfDeterminer["pca"].spatialOrder = 2
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSizeMin = 25

# Final photometry
root.measurement.slots.apFlux = "flux.sinc"
#import lsst.meas.extensions.photometryKron
#import lsst.meas.extensions.rotAngle
#import lsst.meas.extensions.shapeHSM
#root.measurement.algorithms.names += ("rotAngle", "flux.kron",)
#root.measurement.algorithms.names += tuple(["shape.hsm." + s for s in ('bj', 'linear', 'ksb', 'regauss', 
#                                                                       'shapelet')])
#import lsst.meas.extensions.shapeMiyatake
#root.measurement.algorithms.names += ("shape.miyatake.hmfit",)

# Astrometry
root.calibrate.astrometry.solver.filterMap = { 'B': 'g',
                                               'V': 'r',
                                               'R': 'r',
                                               'I': 'i',
                                               'y': 'z',
                                               }

# FH below for QA stuff
root.qa.seeing.fracSrcIni = 0.15
#root.qa.seeing.fwhmIni = 3.465      #pix
#root.qa.seeing.fwhmIni = 4.12      #pix hsc
root.qa.seeing.fwhmMin = 1.5        #pix
root.qa.seeing.fwhmMax = 12.0       #pix
root.qa.seeing.nbinMagHist = 80
root.qa.seeing.magMinHist = -20
root.qa.seeing.magMaxHist = 0.0
root.qa.seeing.nSampleRoughFwhm = 30
####root.qa.seeing.fwhmMarginFinal = 1.5  #pix
root.qa.seeing.fwhmMarginFinal = 0.5  #pix
root.isr.qa.flatness.meshX = 256
root.isr.qa.flatness.meshY = 256
root.isr.qa.flatness.doClip = True
root.isr.qa.flatness.clipSigma = 3.0
root.isr.qa.flatness.nIter = 3
root.isr.qa.doWriteImage.doWriteOssImage = False # "Do we write out overscan-subtracted image?"
root.isr.qa.doWriteImage.doWriteFltImage = False # "Do we write out flatfielded image?"
#root.isr.qa.doWriteImage.doWriteSsbImage = True # "Do we write out flatfielded & sky-subtracted image?"
root.isr.qa.doWriteImage.doDumpSnapshot = True  # "Do we dump snapshot figures of images?"
root.qa.camera = 'suprimecam'
root.isr.qa.camera = 'suprimecam'
root.qa.seeing.camera = 'suprimecam'
