root.calibrate.repair.doCosmicRay = False
#root.calibrate.background.binSize = 1024
root.calibrate.repair.doCrosstalk = False
root.calibrate.repair.doLinearize = False

# Disable background estimation due to ticket #2037
root.calibrate.doBackground = False
root.calibrate.detection.reEstimateBackground = False
root.detection.reEstimateBackground = False


# PSF determination
root.calibrate.detection.thresholdValue = 50
root.calibrate.measurePsf.starSelector.name = "secondMoment"
root.calibrate.measurePsf.psfDeterminer.name = "pca"
root.calibrate.measurePsf.starSelector["secondMoment"].clumpNSigma = 2.0
root.calibrate.measurePsf.psfDeterminer["pca"].nEigenComponents = 4
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSize = 12
root.calibrate.measurePsf.psfDeterminer["pca"].spatialOrder = 4
root.calibrate.measurePsf.psfDeterminer["pca"].kernelSizeMin = 25
root.calibrate.measurePsf.psfDeterminer["pca"].reducedChi2ForPsfCandidates = 10.0
root.calibrate.measurePsf.psfDeterminer["pca"].spatialReject = 10.0


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
