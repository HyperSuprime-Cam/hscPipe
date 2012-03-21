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
#root.doWriteOssImage = True # "Do we write out overscan-subtracted image?"
#root.doWriteFltImage = True # "Do we write out flatfielded image?"
#root.doWriteSsbImage = True # "Do we write out flatfielded & sky-subtracted image?"
#root.doDumpSnapshot = True  # "Do we dump snapshot figures of images?"
    
