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

class SubaruProcessCcdOnsiteTask(SubaruProcessCcdTask):
    """Subclass of SubaruProcessCcdTask, with overrides to turn off heavy
    processing.
    """
    overrides = (applyOnsiteConfigOverrides,)

