#!/usr/bin/env python

import os
import lsst.daf.persistence as dafPersist

def getButler(instrument, rerun=None, **kwargs):
    """Return a butler for the appropriate instrument"""
    if rerun is None:
        rerun = os.getlogin()
    if instrument.lower() in ["hsc"]:
        import lsst.obs.hscSim as obsHsc
        mapper = obsHsc.HscSimMapper(rerun=rerun, **kwargs)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        import lsst.obs.suprimecam as obsSc
        mapper = obsSc.SuprimecamMapper(rerun=rerun, **kwargs)
    else:
        raise RuntimeError("Unrecognised instrument: %s" % instrument)

    return dafPersist.ButlerFactory(mapper=mapper).create()


def getNumCcds(instrument):
    """Return the number of CCDs in an instrument"""
    # XXX This could be done by inspecting the number of Ccds in butler.mapper.camera
    if instrument.lower() in ["hsc"]:
        return 100 # XXX Update when we can handle rotated CCDs
    if instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        return 10
    raise RuntimeError("Unrecognised instrument: %s" % instrument)

