# Overrides for ProcessExposure/ProcessCcd
try:
    processCcd = root.processCcd
except KeyError:
    processCcd = root

# COSMOS catalogue has excess sources around bright stars; increasing numBrightStars cuts through the cruft
processCcd.calibrate.astrometry.solver.numBrightStars = 300

# COSMOS catalogue goes deeper than SDSS
processCcd.calibrate.photocal.magLimit = 24.0
