#!/usr/bin/env python

import os
import lsst.daf.persistence as dafPersist

def getButler(instrument, rerun=os.getlogin()):
    if instrument.lower() in ["hsc"]:
        import lsst.obs.hscSim as obsHsc
        mapper = obsHsc.HscSimMapper(rerun=rerun)
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        import lsst.obs.suprimecam as obsSc
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
    else:
        raise RuntimeError("Unrecognised instrument: %s" % instrument)

    return dafPersist.ButlerFactory(mapper=mapper).create()


def getNumCcds(instrument):
    if instrument.lower() in ["hsc"]:
        return 100 # XXX Update when we can handle rotated CCDs
    if instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        return 10
    raise RuntimeError("Unrecognised instrument: %s" % instrument)

