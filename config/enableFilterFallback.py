# Using this configuration file allows the pipeline to push through
# processing for an image with an unrecognised filter.  While this
# may produce results not suitable for science, it may allow the
# onsite analysis to produce useful results.

try:
    # In case we're being loaded for ProcessExposureTask
    root = root.processCcd
except:
    pass

print "WARNING: loading filter fallbacks; this may produce results not suitable for science"
root.isr.fallbackFilterName = "HSC-I" # Use this for flats
root.calibrate.astrometry.solver.filterMap["UNRECOGNISED"] = "i" # Use this for photocal

# Calibrate directly against i-band: zero color terms
from lsst.meas.photocal.colorterms import ColortermConfig, ColortermGroupConfig
colorterm = ColortermConfig.fromValues("i", "i", 0.0, 0.0, 0.0)
for library in root.calibrate.photocal.colorterms.library.values():
    library.group["UNRECOGNISED"] = colorterm
