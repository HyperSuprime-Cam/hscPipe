root.calibrate.repair.doCosmicRay = False
#root.calibrate.background.binSize = 1024

# Diable background re-estimation due to ticket #2037
root.calibrate.detection.reEstimateBackground = False
root.detection.reEstimateBackground = False
root.calibrate.background.binSize = 500

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
root.measurement.algorithms.names += ("centroid.record",)
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
