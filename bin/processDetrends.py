#!/usr/bin/env python

import os, os.path
import sys
import argparse
import mpi4py.MPI as mpi
import pbasf2 as pbasf

import lsst.daf.base as dafBase
import lsst.afw.image as afwImage
from lsst.pipe.base import Struct
import hsc.pipe.base.camera as hscCamera
import hsc.pipe.tasks.processCcd as hscProcessCcd
import hsc.pipe.tasks.detrends as hscDetrends



def processDetrends(instrument, rerun, detrend, frameList, outRoot=None, outName=None):
    comm = mpi.COMM_WORLD

    if outName is None:
        outName = detrend + "-" + rerun + "-%d.fits"

    butler = hscCamera.getButler(instrument, rerun)
    numCcds = hscCamera.getNumCcds(instrument)
    ccdList = range(numCcds)# if False else range(2)

    config = hscDetrends.DetrendConfig()
    override = os.path.join(os.environ['HSCPIPE_DIR'], 'config', "%s-%s.py" % (instrument, detrend))
    try:
        config.load(override)
    except Exception, e:
        sys.stderr.write("Unable to read override file (%s): %s\n" % (override, e))

    # Scatter/gather: process individual detrends, get backgrounds
    worker = Worker(butler, config, detrend, frameList, outName)
    backgrounds = pbasf.ScatterJob(comm, worker.process, ccdList, root=0)

    if detrend.lower() != "flat":
        # All done!
        return

    # Together: calculate scalings
    workList = pbasf.SafeCall(worker.scales, backgrounds) if comm.Get_rank() == 0 else None

    # Scatter result: combine
    pbasf.ScatterJob(comm, worker.combineFlat, workList, root=0)


class Worker(object):
    def __init__(self, butler, config, detrend, frameList, outName, version="000"):
        self.butler = butler
        self.task = hscDetrends.DetrendTask(config=config)
        self.detrend = detrend.lower()
        self.frameList = frameList
        self.outName = outName
        self.version = version
        self.maskData = []

    def _getDataRef(self, frame, ccd):
        dataId = {'visit': frame, 'ccd': ccd}
        dataRefList = [ref for ref in self.butler.subset(datasetType='raw', **dataId)]
        if len(dataRefList) != 1:
            raise RuntimeError("Unable to get data reference for dataId=%s" % dataId)
        return dataRefList[0]

    def process(self, ccd):
        results = []
        for frame in self.frameList:
            try:
                result = self.processSingle(frame, ccd)
            except Exception, e:
                sys.stderr.write("Failed to process frame %d ccd %d: %s\n" % (frame, ccd, e))
                continue
            if self.detrend == "mask":
                results.append(result)
            elif self.detrend == "flat":
                results.append(result.background)

        if self.detrend == "flat":
            # CCDs are inter-dependent: need to ensure consistent scaling
            return results

        # All CCDs are independent, so we can combine now
        if self.detrend == "mask":
            self.mask(ccd, results)
        else:
            self.combine(ccd)

    def processSingle(self, frame, ccd):
        dataRef = self._getDataRef(frame, ccd)
        result = self.task.process.run(self.detrend, dataRef)

        if self.detrend != "mask":
            # Save exposure for combination
            dataRef.put(result.exposure, "calexp")
            del result.exposure

        return result

    def scales(self, backgrounds):
        assert self.detrend == "flat"
        import numpy
        result = self.task.scale.run(numpy.array(backgrounds, dtype=float))
        ccdScales = result.components
        expScales = result.exposures
        return [Struct(ccd=ccd, ccdScale=scale, expScales=expScales) for ccd, scale in enumerate(ccdScales)]

    def combineFlat(self, work):
        self.combine(work.ccd, expScales=work.expScales, ccdScale=work.ccdScale)

    def combine(self, ccd, expScales=None, ccdScale=None):
        dataRefList = []
        for frame in self.frameList:
            try:
                dataRefList.append(self._getDataRef(frame, ccd))
            except:
                pass
        combined = self.task.combine.run(dataRefList, expScales=expScales, finalScale=ccdScale)
        self.write(ccd, combined)

    def mask(self, ccd, processResults):
        footprints = [data.footprintSets for data in processResults]
        dimensions = [data.dim for data in processResults]
        combined = self.task.mask.run(footprints, dimensions)
        self.write(ccd, combined)

    def write(self, ccd, image):
        filterName = None
        midTime = 0
        for frame in self.frameList:
            dataId = {'visit': frame, 'ccd': ccd}
            md = self.butler.get("calexp_md", dataId)

            calib = afwImage.Calib(md)
            midTime += calib.getMidTime().mjd()

            if self.detrend in ("flat", "fringe"):
                thisFilter = self.butler.queryMetadata("raw", "filter", format="filter", dataId=dataId)[0]
                if filterName is None:
                    filterName = thisFilter
                elif filterName != thisFilter:
                    print "WARNING: Filter mismatch for %s: %s vs %s" % (dataId, thisFilter, filterName)

        midTime /= len(self.frameList)
        date = str(dafBase.DateTime(midTime, dafBase.DateTime.MJD).toPython().date())

        self.butler.put(image, self.detrend, filter=filterName, mystery=self.version, num=0,
                        calibDate=date, ccd=ccd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--instrument", type=str, required=True, help="Instrument name")
    parser.add_argument("-r", "--rerun", type=str, required=True, help="Rerun name")
    parser.add_argument("-d", "--detrend", type=str, choices=["bias", "dark", "flat", "fringe", "mask"],
                        help="Detrend type", required=True)
    parser.add_argument("-o", "--out", type=str, help="Pattern for output name, with %%d for ccd number")
    parser.add_argument("frames", type=int, nargs="+", help="Frames to combine")

    args = parser.parse_args()

    processDetrends(args.instrument, args.rerun, args.detrend, args.frames, args.out)
