from lsst.obs.subaru.isr import SuprimeCamIsrTask
root.process.isr.retarget(SuprimeCamIsrTask)
root.process.isr.doBias = False
root.process.isr.doDark = False
root.process.isr.doFlat = True

# Mask out bright sources
root.process.doDetection = True
root.process.detection.thresholdValue = 5.0
root.process.detection.includeThresholdMultiplier = 5.0
root.combine.maskDetected = True
