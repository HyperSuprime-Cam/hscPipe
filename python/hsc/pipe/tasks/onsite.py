#!/usr/bin/env python

import lsst.pex.config as pexConfig

from .processCcd import SubaruProcessCcdTask
from .qa import QaTask

def applyOnsiteConfigOverrides(root):
    root.isr.qa.doWriteOss = True
    root.isr.qa.doThumbnailOss = True
    root.isr.qa.doWriteFlattened = True
    root.isr.qa.doThumbnailFlattened = True

class SubaruProcessCcdOnsiteConfig(SubaruProcessCcdTask.ConfigClass):
    qa = pexConfig.ConfigurableField(target = QaTask, doc = "QA analysis")

class SubaruProcessCcdOnsiteTask(SubaruProcessCcdTask):
    """Subclass of SubaruProcessCcdTask with additional quick-look QA analysis
    (note that some is already present in lower-level packages, especially the
    Subaru ISR tasks in obs_subaru).
    """
    ConfigClass = SubaruProcessCcdOnsiteConfig

    overrides = (applyOnsiteConfigOverrides,)

    def __init__(self, *args, **kwargs):
        SubaruProcessCcdTask.__init__(self, *args, **kwargs)
        self.makeSubtask("qa")

    def run(self, sensorRef):
        result = SubaruProcessCcdTask.run(self, sensorRef)
        self.qa.run(sensorRef, result.exposure, result.sources)
        return result
