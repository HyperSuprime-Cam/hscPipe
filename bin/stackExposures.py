#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal
import time
import optparse

import lsst.pipette.runHsc as runHsc
import lsst.obs.hscSim as obsHsc
import lsst.obs.suprimecam as obsSc
import lsst.pipette.readwrite as pipReadWrite

import hsc.meas.mosaic.mosaic as hscMosaic
import hsc.meas.mosaic.stack as hscStack

def sigalrm_handler(signum, frame):
    sys.stderr.write('Signal handler called with signal %s\n' % (signum))
signal.signal(signal.SIGALRM, sigalrm_handler)

def main():
    parser = optparse.OptionParser()
    parser.add_option("-r", "--rerun",
                      type=str, default=None,
                      help="rerun name to take corrected frames from and write stack images to.")
    parser.add_option("-i", "--instrument",
                      type=str, default='hsc',
                      help="instument to treat (hsc or suprimecam)")
    parser.add_option("-p", "--program",
                      type=str, default=None,
                      help="program name (e.g. COSMOS_0)")
    parser.add_option("-f", "--filter",
                      type=str, default=None,
                      help="filter name (e.g. W-S-I+)")
    parser.add_option("-d", "--dateObs",
                      type=str, default=None,
                      help="(optional) dataObs (e.g. 2008-11-27)")
    parser.add_option("-w", "--workDirRoot",
                      type=str, default=".",
                      help="root working directory (working dir will be root/program/filter)")
    parser.add_option("-s", "--destWcs",
                      type=str, default=None,
                      help="destination wcs")
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args

    print "rerun=%s, instrument=%s, program=%s, filter=%s, dateObs=%s, workDirRoot=%s, destWcs=%s, args=%s " % \
        (opts.rerun, opts.instrument, opts.program, opts.filter, opts.dateObs, opts.workDirRoot, opts.destWcs, sys.argv)

    try:
        ProcessMosaicStack(rerun=opts.rerun, instrument=opts.instrument, program=opts.program, filter=opts.filter, dateObs=opts.dateObs, workDirRoot=opts.workDirRoot, destWcs=opts.destWcs)
        return 0;
    except:
        pbasf.ReportError("Total catastrophic failure")
        print "THIS ERROR SHALL NOT HAVE APPEARED."
        mpi.COMM_WORLD.Abort(1)
        return 1
        
def ProcessMosaicStack(rerun=None, instrument=None, program=None, filter=None, dateObs=None, workDirRoot=None, destWcs=None):
    if instrument.lower() in ["hsc"]:
        mapper = obsHsc.HscSimMapper(rerun=rerun)
        nCCD = 100
    elif instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
        mapper = obsSc.SuprimecamMapper(rerun=rerun)
        nCCD = 10
    else:
        raise RuntimeError("unknown instrument: %s" % (instrument))

    ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})

    if (dateObs == None):
        lFrameId = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter))
        lPointing = ioMgr.inButler.queryMetadata('calexp', None, 'pointing', dict(field=program, filter=filter))
    else:
        lFrameId = ioMgr.inButler.queryMetadata('calexp', None, 'visit', dict(field=program, filter=filter, dateObs=dateObs))
        lPointing = ioMgr.inButler.queryMetadata('calexp', None, 'pointing', dict(field=program, filter=filter, dateObs=dateObs))
    print lFrameId
    print lPointing

    config = {"filter":filter, \
              "stackId":lPointing[0], \
              "program":program, \
              "dateObs":dateObs, \
              "subImgSize":4096, \
              "imgMargin":256, \
              "fileIO":True, \
              "writePBSScript":False, \
              "skipMosaic":False, \
              "workDirRoot":workDirRoot}

    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()

    # create ccdId's
    lCcdId = range(nCCD)

    indexes = []
    if rank == 0:
        # phase 1
        lFrameIdExist = pbasf.SafeCall(phase1, ioMgr, lFrameId, lCcdId, workDirRoot)

        # phase 2
        nx, ny = pbasf.SafeCall(phase2, ioMgr, lFrameIdExist, lCcdId, instrument, rerun, destWcs, config)

        print 'nx = ', nx, ' ny = ', ny

        indexes = [(ix, iy) for ix in range(nx) for iy in range(ny)]
        comm.bcast(indexes, root=0)
    else:
        indexes = comm.bcast(indexes, root=0)

    # phase 3
    if rank == 0:
        phase3 = None
    else:
        phase3 = Phase3Worker(rerun=rerun, instrument=instrument, config=config)
    pbasf.ScatterJob(comm, phase3, [index for index in indexes], root=0)

    if rank == 0:
        # phase 4
        pbasf.SafeCall(phase4, ioMgr, instrument, rerun, config)

def phase1(ioMgr, lFrameId, lCcdId, workDirRoot):
    if True:
        return hscMosaic.mosaic(ioMgr, lFrameId, lCcdId, outputDir=workDirRoot)
    else:
        lFrameIdExist = []
        for frameId in lFrameId:
            good = True
            for ccdId in lCcdId:
                good |= ioMgr.inButler.datasetExists('calexp', dict(visit=frameId, ccd=ccdId))
                good |= ioMgr.inButler.datasetExists('wcs', dict(visit=frameId, ccd=ccdId))
                good |= ioMgr.inButler.datasetExists('fcr', dict(visit=frameId, ccd=ccdId))
                if not good:
                    break
            if good:
                lFrameIdExist.append(frameId)
        return lFrameIdExist

def phase2(ioMgr, lFrameId, lCcdId, instrument, rerun, destWcs, config):
    fileList = []
    for frameId in lFrameId:
        for ccdId in lCcdId:
            try:
                fname = ioMgr.read('calexp_filename', dict(visit=frameId, ccd=ccdId))[0][0]
            except Exception, e:
                print "failed to get file for %s:%s" % (frameId, ccdId)
                continue
            if os.path.isfile(fname):
                fileList.append(fname)
            else:
                print "file %s does not exist " % (fname)

    subImgSize = config['subImgSize']
    imgMargin = config['imgMargin']
    fileIO = config['fileIO']
    writePBSScript = config['writePBSScript']
    skipMosaic = config['skipMosaic']
    program = config['program']
    filter = config['filter']
    dateObs = config['dateObs']
    workDirRoot = config['workDirRoot']

    workDir = os.path.join(workDirRoot, program, filter)
    try:
        os.makedirs(workDir)
    except OSError:
        print "Working directory already exists"

    if destWcs != None:
        destWcs = os.path.abspath(destWcs)

    return hscStack.stackInit(ioMgr,
                              fileList, subImgSize, imgMargin,
                              fileIO, writePBSScript,
                              workDir=workDir, skipMosaic=skipMosaic,
                              rerun=rerun, instrument=instrument,
                              program=program, filter=filter, dateObs=dateObs,
                              destWcs=destWcs)

class Phase3Worker:
    def __init__(self, rerun=None, instrument="hsc", config=None):
        self.rerun = rerun
        self.instrument = instrument
        self.config = config
    
    def __call__(self, t_ix_iy):
        if self.instrument.lower() in ["hsc"]:
            mapper = obsHsc.HscSimMapper(rerun=self.rerun)
        elif self.instrument.lower() in ["suprimecam", "suprime-cam", "sc"]:
            mapper = obsSc.SuprimecamMapper(rerun=self.rerun)

        ioMgr = pipReadWrite.ReadWrite(mapper, ['visit', 'ccd'], config={})

        ix = t_ix_iy[0]
        iy = t_ix_iy[1]
        print "Started processing %d,%d in %s, %d" % (ix, iy, os.uname()[1], os.getpid())

        stackId = self.config['stackId']
        program = self.config['program']
        filter = self.config['filter']
        dateObs = self.config['dateObs']
        subImgSize = self.config['subImgSize']
        imgMargin = self.config['imgMargin']
        fileIO = self.config['fileIO']
        skipMosaic = self.config['skipMosaic']
        workDirRoot = self.config['workDirRoot']
        workDir = os.path.join(workDirRoot, program, filter)

        hscStack.stackExec(ioMgr, ix, iy, stackId, subImgSize, imgMargin, fileIO=fileIO, workDir=workDir, skipMosaic=skipMosaic, filter=filter)

def phase4(ioMgr, instrument, rerun, config):
    stackId = config['stackId']
    program = config['program']
    filter = config['filter']
    subImgSize = config['subImgSize']
    imgMargin = config['imgMargin']
    fileIO = config['fileIO']
    workDirRoot = config['workDirRoot']
    workDir = os.path.join(workDirRoot, program, filter)

    hscStack.stackEnd(ioMgr, stackId, subImgSize, imgMargin, fileIO=fileIO,
                      workDir=workDir, filter=filter)

if __name__ == "__main__":
    print "argv=", sys.argv
    main()
