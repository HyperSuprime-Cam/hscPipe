from lsst.obs.subaru.isr import SuprimeCamIsrTask
root.process.isr.retarget(SuprimeCamIsrTask)
root.process.isr.doBias = False
root.process.isr.doDark = False
root.process.isr.doFlat = False
root.process.isr.doCrosstalk = False
root.process.isr.doGuider = False

root.process.doDetection = False
root.combine.maskDetected = False
