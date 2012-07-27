root.postprocess.doQa = True
root.isr.qa.doWriteOss = True
root.isr.qa.doThumbnailOss = True
root.isr.qa.doWriteFlattened = True
root.isr.qa.doThumbnailFlattened = True

# Disable model mags becaause they are slow (could disable more than this if desired).
root.measurement.algorithms.names = [name for name in root.measurement.algorithms.names
                                     if not name.startswith("multishapelet")]
root.measurement.slots.modelFlux = "flux.gaussian"
