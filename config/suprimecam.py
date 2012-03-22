root.doWriteIsr = True
##root.doWriteIsr = False
root.isr.methodList=["doConversionForIsr", "doSaturationDetection",
#                     "doOverscanCorrection", "doVariance", "doFlatCorrection"]
                     "doOverscanCorrectionQa", "doVariance", "doWriteOssImageQa", "doFlatCorrectionQa", "doWriteFltImageQa"]
##root.isr.doWrite = False
##root.isr.doWrite = True

root.calibrate.repair.doCosmicRay = True
root.calibrate.repair.cosmicray.nCrPixelMax = 1000000
root.calibrate.background.binSize = 1024

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

# Astrometry
root.calibrate.astrometry.solver.filterMap = { 'B': 'g',
                                               'V': 'r',
                                               'R': 'r',
                                               'I': 'i',
                                               'y': 'z',
                                               }

# FH below for QA stuff
root.qa.seeing.fwhmIni = 3.465
root.qa.seeing.fwhmMin = 1.5
root.qa.seeing.fwhmMax = 12.0
root.qa.seeing.nbinMagHist = 80
root.qa.seeing.magMinHist = -20
root.qa.seeing.magMaxHist = 0.0
root.qa.seeing.nSampleRoughFwhm = 30
root.qa.seeing.fwhmMarginFinal = 1.5
root.qa.flatness.meshX = 256
root.qa.flatness.meshY = 256
root.qa.flatness.doClip = True
root.qa.flatness.clipSigma = 3.0
root.qa.flatness.nIter = 3
root.qa.doWriteOssImage = True # "Do we write out overscan-subtracted image?"
root.qa.doWriteFltImage = True # "Do we write out flatfielded image?"
#root.qa.doWriteSsbImage = True # "Do we write out flatfielded & sky-subtracted image?"
root.qa.doDumpSnapshot = True  # "Do we dump snapshot figures of images?"
    
