#!/usr/bin/env python

def write(self, butler, dataId, struct, wcs=None):
    if wcs is None:
        wcs = struct.exposure.getWcs()
        self.log.log(self.log.WARN, "WARNING: No new WCS provided")

    # Apply WCS to sources
    # No longer handling matchSources explicitly - these should all be in calib.sources,
    # or there's a bug in the calibrate task.
    struct.exposure.setWcs(wcs)
    for sources in (struct.sources, struct.calib.sources):
        if sources is None:
            continue
        for s in sources:
            s.updateCoord(wcs)

    self.writeMatches(struct.calib.matches, struct.calib.matchMeta)

    butler.put(struct.exposure, 'calexp', dataId)
    butler.put(struct.sources, 'src', dataId)
    butler.put(normalizedMatches, 'icMatch', dataId)
    butler.put(struct.calib.psf, 'psf', dataId)
    butler.put(struct.calib.apCorr, 'apCorr', dataId)
    butler.put(struct.calib.sources, 'icSrc', dataId)
