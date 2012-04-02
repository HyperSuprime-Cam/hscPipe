#!/usr/bin/env python

import sys
import argparse
import mpi4py.MPI as mpi
import pbasf2 as pbasf

import hsc.pipe.base.camera as hscCamera
import hsc.pipe.tasks.processCcd as hscProcessCcd
import hsc.pipe.tasks.detrends as hscDetrends



def processDetrends(instrument, rerun, detrend, frameList, outName=None):
    comm = mpi.COMM_WORLD

    if outName is None:
        outName = "detrend-%s-\%d.fits" % rerun

    butler = hscCamera.getButler(instrument, rerun)
    numCcds = hscCamera.getNumCcds(instrument)
    ccdList = range(numCcds)

    # Scatter/gather: process individual detrends, get backgrounds
    worker = Worker(butler, config, detrend, frameList, outName)
    results = pbasf.ScatterJob(comm, worker.run, ccdList, root=0)

    # Together: 
    if comm.Get_rank() == 0:
        backgrounds = pbasf.SafeCall(worker.scales, results)

    # Scatter result: combine
    pbasf.QueryToRoot(comm, worker.combine, lambda ccd: results[ccd], ccdList, root=0)


class Worker(object):
    def __init__(self, butler, config, detrend, frameList, outName):
        self.butler = butler
        self.task = hscDetrends.DetrendTask(config=config)
        self.detrend = detrend.lower()
        self.frameList = frameList
        self.outName = outName
        self.maskData = []

    def process(self, ccd):
        backgrounds = []
        for frame in frameList:
            dataId = {'visit': frame, 'ccd': ccd}
            dataRefList = [ref for ref in self.butler.subset(datasetType='raw', **dataId)]
            assert len(dataRefList) == 1
            dataRef = dataRefList[0]

            try:
                result = self.task.process.run(self.detrend, dataRef)
            except Exception, e:
                sys.stderr.write("Failed to process %s: %s\n" % (dataId, e))
                raise

            if self.detrend == "mask":
                self.maskData.append(result)
            else:
                self.butler.put("calexp", dataId, result.exposure)

            backgrounds.append(result.background)

        return backgrounds

    def scales(self, backgrounds):
        if detrend == "flat":
            import numpy
            ccdScales, expScales = self.task.scale.run(numpy.array(backgrounds))
        else:
            ccdScales = [None] * len(backgrounds)
            expScales = None
        return ccdScales, expScales

    def combine(self, ccd, scales):
        expScales = scales.expScales
        ccdScale = scales.ccdScale

        if self.detrend == "mask":
            footprints = [data.footprintSets for data in self.maskData]
            dimensions = [data.dim for data in self.maskData]
            combined = self.task.mask.run(footprints, dimensions)
        else:
            combined = self.task.combine.run(dataRefList, expScales=expScales, finalScale=ccdScale)
            
        combined.writeFits(outName % ccd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--instrument", type=str, required=True, help="Instrument name")
    parser.add_argument("-r", "--rerun", type=str, required=True, help="Rerun name")
    parser.add_argument("-d", "--detrend", type=str, choices=["bias", "dark", "flat", "fringe", "mask"],
                        help="Detrend type")
    parser.add_argument("-o", "--out", type=str, help="Pattern for output name, with \%d for ccd number")
    parser.add_argument("frames", type=int, nargs="+", help="Frames to combine")

    args = parser.parse_args()

    processDetrends(args.instrument, args.rerun, args.detrend, args.frames, args.out)
