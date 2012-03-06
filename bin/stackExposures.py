#!/usr/bin/env python

import sys
import mpi4py.MPI as mpi
import pbasf2 as pbasf
import os
import signal
import time
import optparse

import lsst.pex.config as pexConfig
import lsst.obs.hscSim as obsHsc
import lsst.obs.suprimecam as obsSc

import hsc.pipe.base.camera as hscCamera
import hsc.meas.mosaic.mosaic as hscMosaic
import hsc.meas.mosaic.config as hscMC
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
    parser.add_option("-m", "--doMatchPsf",
		      default=False, action='store_true',
		      help="match PSFs before stacking (default=%default)")
    
    (opts, args) = parser.parse_args()

    if not opts.rerun or not opts.program or not opts.filter:
        parser.print_help()
        raise SystemExit("failed to parse arguments")

    sys.argv = [sys.argv[0]] + args

    print "rerun=%s, instrument=%s, program=%s, filter=%s, dateObs=%s, workDirRoot=%s, destWcs=%s, doMatchPsf=%s  args=%s " % \
        (opts.rerun, opts.instrument, opts.program, opts.filter, opts.dateObs, opts.workDirRoot, opts.destWcs, str(opts.doMatchPsf), sys.argv)

    try:
        ProcessMosaicStack(rerun=opts.rerun, instrument=opts.instrument, program=opts.program, filter=opts.filter, dateObs=opts.dateObs, workDirRoot=opts.workDirRoot, destWcs=opts.destWcs, doMatchPsf=opts.doMatchPsf)
        return 0;
    except:
        pbasf.ReportError("Total catastrophic failure")
        print "THIS ERROR SHALL NOT HAVE APPEARED."
        mpi.COMM_WORLD.Abort(1)
        return 1
        
def ProcessMosaicStack(rerun=None, instrument=None, program=None, filter=None,
                       dateObs=None, workDirRoot=None, destWcs=None, doMatchPsf=False):
    butler = hscCamera.getButler(instrument, rerun)
    nCCD = hscCamera.getNumCcds(instrument)

    dataId = dict(field=program, filter=filter)
    if dateObs is not None:
        dataId['dateObs'] = dateObs

    lFrameId = butler.queryMetadata('calexp', None, 'visit', dataId)
    lPointing = butler.queryMetadata('calexp', None, 'pointing', dataId)

    print lFrameId
    print lPointing

    mosaicConfig = hscMC.HscMosaicConfig()
    stackConfig = hscMC.HscStackConfig()
    stackConfig.filterName = filter
    stackConfig.stackId = lPointing[0]
    stackConfig.program = program
    stackConfig.dateObs = dateObs
    stackConfig.workDirRoot = workDirRoot

    comm = mpi.COMM_WORLD
    rank = comm.Get_rank()

    # create ccdId's
    lCcdId = range(nCCD)

    dataPack = {
        'indexes' : [],
        'fileList' : [],
        'wcs' : None
        }
    #indexes = []
    if rank == 0:
        # phase 1
        lFrameIdExist = pbasf.SafeCall(phase1, butler, lFrameId, lCcdId, workDirRoot, mosaicConfig)
        print lFrameIdExist

        # phase 2
        nx, ny, fileList, wcs = pbasf.SafeCall(phase2, butler, lFrameIdExist, lCcdId, instrument, rerun,
                                               destWcs, stackConfig)

        print 'nx = ', nx, ' ny = ', ny

        indexes = [(ix, iy) for ix in range(nx) for iy in range(ny)]
        dataPack['indexes'] = indexes
        dataPack['fileList'] = fileList
        dataPack['wcs'] = wcs
        
        comm.bcast(dataPack, root=0)
    else:
        dataPack = comm.bcast(dataPack, root=0)
        indexes = dataPack['indexes']
        fileList = dataPack['fileList']
        wcs = dataPack['wcs']

    # phase 3 (measure PSFs in warped images)
    sigmas = []
    if doMatchPsf:
        if rank == 0:
            phase3a = None
        else:
            phase3a = Phase3aWorker(butler, config=stackConfig, wcs=wcs)
        sigmas = pbasf.ScatterJob(comm, phase3a, [f for f in fileList], root=0)

    # phase 3b
    dummy = None
    if rank == 0: # or sigmas is None:
        phase3b = None
        comm.bcast(sigmas, root=0)
    else:
        sigmas = comm.bcast(sigmas, root=0)
        matchPsf = None
        print "rank/sigmas:", rank, type(sigmas), sigmas
        if sigmas:
            maxSigma = max(sigmas)
            sigma1 = maxSigma
            sigma2 = 2.0*maxSigma
            kwid = int(4.0*sigma2) + 1
            peakRatio = 0.1
            matchPsf = ['DoubleGaussian', kwid, kwid, sigma1, sigma2, peakRatio]

        phase3b = Phase3bWorker(butler, config=stackConfig, matchPsf=matchPsf)
    pbasf.ScatterJob(comm, phase3b, [index for index in indexes], root=0)

    
    if rank == 0:
        # phase 4
        pbasf.SafeCall(phase4, butler, instrument, rerun, stackConfig)

def phase1(butler, lFrameId, lCcdId, workDirRoot, mosaicConfig):
    if True:
        return hscMosaic.mosaic(butler, lFrameId, lCcdId, mosaicConfig, outputDir=workDirRoot)
    else:
        lFrameIdExist = []
        for frameId in lFrameId:
            good = True
            for ccdId in lCcdId:
                good |= butler.datasetExists('calexp', dict(visit=frameId, ccd=ccdId))
                good |= butler.datasetExists('wcs', dict(visit=frameId, ccd=ccdId))
                good |= butler.datasetExists('fcr', dict(visit=frameId, ccd=ccdId))
                if not good:
                    break
            if good:
                lFrameIdExist.append(frameId)
        return lFrameIdExist

def phase2(butler, lFrameId, lCcdId, instrument, rerun, destWcs, config):
    fileList = []
    for frameId in lFrameId:
        for ccdId in lCcdId:
            try:
                fname = butler.get('calexp_filename', dict(visit=frameId, ccd=ccdId))[0]
            except Exception, e:
                print "failed to get file for %s:%s" % (frameId, ccdId)
                continue
            if os.path.isfile(fname):
                fileList.append(fname)
            else:
                print "file %s does not exist " % (fname)

    workDir = os.path.join(config.workDirRoot, config.program, config.filterName)
    try:
        os.makedirs(workDir)
    except OSError:
        print "Working directory already exists"

    if destWcs != None:
        destWcs = os.path.abspath(destWcs)

    return hscStack.stackInit(butler,
                              fileList, config.subImageSize, config.imageMargin,
                              config.fileIO, config.writePbsScript,
                              workDir=workDir, skipMosaic=config.skipMosaic,
                              rerun=rerun, instrument=instrument,
                              program=config.program, filter=config.filterName,
                              dateObs=config.dateObs, destWcs=destWcs)

class Phase3Worker:
    def __init__(self, butler, config, wcs=None):
        self.butler = butler
        self.config = config
        self.wcs = wcs
        
    def __call__(self, fname):
        print "Started measuring warped PSF for %s in %s, %d" % (fname, os.uname()[1], os.getpid())

        stackId = self.config['stackId']
        program = self.config['program']
        filter = self.config['filter']
        dateObs = self.config['dateObs']
        subImgSize = self.config['subImgSize']
        fileIO = self.config['fileIO']
        skipMosaic = self.config['skipMosaic']
        workDirRoot = self.config['workDirRoot']
        workDir = os.path.join(workDirRoot, program, filter)

        return hscStack.stackMeasureWarpedPsf(fname, self.wcs, butler=butler, fileIO=True,
                                              skipMosaic=skipMosaic)



class Phase3bWorker:
    def __init__(self, butler, config=None, matchPsf=None):
        self.butler = butler
        self.config = config
        self.matchPsf = matchPsf
    
    def __call__(self, t_ix_iy):
        ix = t_ix_iy[0]
        iy = t_ix_iy[1]
        print "Started processing %d,%d in %s, %d" % (ix, iy, os.uname()[1], os.getpid())

        stackId = self.config.stackId
        program = self.config.program
        filter = self.config.filterName
        dateObs = self.config.dateObs
        subImgSize = self.config.subImageSize
        imgMargin = self.config.imageMargin
        fileIO = self.config.fileIO
        skipMosaic = self.config.skipMosaic
        workDirRoot = self.config.workDirRoot
        workDir = os.path.join(workDirRoot, program, filter)

        hscStack.stackExec(self.butler, ix, iy, stackId, subImgSize, imgMargin, fileIO=fileIO,
                           workDir=workDir, skipMosaic=skipMosaic, filter=filter, matchPsf=self.matchPsf)

def phase4(butler, config):
    stackId = config.stackId
    program = config.program
    filter = config.filterName
    subImgSize = config.subImageSize
    imgMargin = config.imageMargin
    fileIO = config.fileIO
    workDirRoot = config.workDirRoot
    workDir = os.path.join(workDirRoot, program, filter)

    hscStack.stackEnd(butler, stackId, subImgSize, imgMargin, fileIO=fileIO,
                      workDir=workDir, filter=filter)

if __name__ == "__main__":
    print "argv=", sys.argv
    main()
