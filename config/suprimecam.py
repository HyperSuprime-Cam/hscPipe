root.calibrate.repair.doCosmicRay = True
root.calibrate.repair.cosmicray.nCrPixelMax = 1000000
root.calibrate.repair.cosmicray.cond3_fac2 = 0.4
root.calibrate.background.binSize = 1024
root.calibrate.detection.background.binSize = 1024

# PSF determination
##root.calibrate.measurePsf.starSelector.name = "secondMoment"
import lsst.meas.astrom.catalogStarSelector
root.calibrate.measurePsf.starSelector.name = "catalog"
root.calibrate.measurePsf.psfDeterminer.name = "pca"
root.calibrate.measurePsf.starSelector["secondMoment"].fluxLim = 60000.0
root.calibrate.measurePsf.starSelector["secondMoment"].clumpNSigma = 2.0
root.calibrate.measurePsf.psfDeterminer["pca"].nEigenComponents = 6
root.calibrate.measurePsf.psfDeterminer["pca"].kernelScaling = 10
root.calibrate.measurePsf.psfDeterminer["pca"].spatialOrder = 4
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSizeMin = 25
root.calibrate.measurePsf.psfDeterminer["pca"].constantWeight = False

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
