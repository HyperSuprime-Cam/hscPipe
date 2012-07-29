#!/usr/bin/env python

import lsst.pex.config as pexConfig

from .processCcd import SubaruProcessCcdTask
from .qa import QaTask

def applyOnsiteConfigOverrides(root):
    root.isr.qa.doWriteOss = True
    root.isr.qa.doThumbnailOss = True
    root.isr.qa.doWriteFlattened = True
    root.isr.qa.doThumbnailFlattened = True
    # Disable these just because they're slow and we don't need them for quick-look.
    # Could probably remove even more, but I don't know what the quick-look requires.
    root.measurement.algorithms.names.discard("multishapelet.psf")
    root.measurement.algorithms.names.discard("multishapelet.exp")
    root.measurement.algorithms.names.discard("multishapelet.dev")
    root.measurement.algorithms.names.discard("multishapelet.combo")
    root.measurement.algorithms.names.discard("centroid.gaussian") 
    root.measurement.slots.modelFlux = "flux.gaussian"

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
