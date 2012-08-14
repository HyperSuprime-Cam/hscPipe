# Disable deblender
root.measurement.doReplaceWithNoise = False
root.doDeblend = False

# Disable multiShapelet
try:
    import lsst.meas.extensions.multiShapelet
    for alg in lsst.meas.extensions.multiShapelet.algorithms:
        root.measurement.algorithms.names.discard(alg)
    root.measurement.slots.modelFlux = "flux.gaussian"
except ImportError:
    pass
