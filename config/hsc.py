root.doWriteIsr = False
root.isr.methodList=["doConversionForIsr", "doSaturationDetection",
                     "doOverscanCorrection", "doVariance", "doFlatCorrection"]
root.isr.doWrite = False

root.calibrate.repair.doCosmicRay = True
root.calibrate.repair.cosmicray.nCrPixelMax = 100000
root.calibrate.background.binSize = 1024

# PSF determination
root.calibrate.measurePsf.starSelector.name = "secondMoment"
root.calibrate.measurePsf.psfDeterminer.name = "pca"
root.calibrate.measurePsf.starSelector["secondMoment"].clumpNSigma = 2.0
root.calibrate.measurePsf.psfDeterminer["pca"].nEigenComponents = 4
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSize = 7
root.calibrate.measurePsf.psfDeterminer["pca"].spatialOrder = 2
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSizeMin = 25

# Final photometry
root.measurement.slots.apFlux = "flux.sinc"

# Astrometry
import hsc.pipe.tasks.distortion
root.calibrate.astrometry.distortion.name = "hsc"
