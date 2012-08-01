from lsst.obs.subaru.isr import SuprimeCamIsrTask
root.process.isr.retarget(SuprimeCamIsrTask)
root.process.isr.doBias = False
root.process.isr.doDark = False
root.process.isr.doFlat = True

root.process.detection.thresholdValue = 3.0
root.process.detection.includeThresholdMultiplier = 1.0

root.mask.maskFraction = 0.5
