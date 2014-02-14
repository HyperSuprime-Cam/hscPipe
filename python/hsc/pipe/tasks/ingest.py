from lsst.pex.config import Config, ConfigurableField
from lsst.obs.subaru.ingest import HscIngestTask
from hsc.pipe.base.pool import Pool, abortOnError
from hsc.pipe.base.pbs import PbsPoolTask

class PoolIngestTask(PbsPoolTask, HscIngestTask):
    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numCcds = sum(1 for raft in parsedCmd.butler.get("camera") for ccd in afwCg.cast_Raft(raft))
        numCycles = int(math.ceil(numCcds/float(numNodes*numProcs)))
        numExps = len(cls.RunnerClass.getTargetList(parsedCmd))
        return time*numExps*numCycles

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doPbs = kwargs.pop("doPbs", False)
        parser = ArgumentParser(name="processExposure", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", level="visit",
                               help="data ID, e.g. --id visit=12345")
        return parser

    @abortOnError
    def run(self, args):
        pool = Pool(None)

        # Read and move files in parallel
        filenameList = args.files
        del args.files # Can be large, and nodes don't need to know
        infoList = pool.map(self.runFile, filenameList, args)

        # Stuff database in serial
        context = self.register.openRegistry(args.butler, create=args.create, dryrun=args.dryrun)
        with context as registry:
            for hduInfoList in infoList:
                if hduInfoList is None:
                    continue
                for info in hduInfoList:
                    self.register.addRow(registry, info, dryrun=args.dryrun, create=args.create)
            self.register.addVisits(registry, dryrun=args.dryrun)

    def runFile(self, infile, args):
        """Examine and ingest a single file

        This method only runs on the slave nodes.

        @param infile: File to process
        @param args: Parsed command-line arguments
        @return parsed information from FITS HDUs
        """
        if self.isBadFile(infile, args.badFile):
            self.log.warn("Skipping declared bad file %s" % infile)
            return None
        fileInfo, hduInfoList = self.parse.getInfo(infile)
        if self.isBadId(fileInfo, args.badId.idList):
            self.log.warn("Skipping declared bad file %s: %s" % (infile, fileInfo))
            return None
        outfile = self.parse.getDestination(args.butler, fileInfo, infile)
        # Note: not doing "self.register.check(registry, fileInfo)" --- it's serial, and makes no difference
        ingested = self.ingest(infile, outfile, mode=args.mode, dryrun=args.dryrun)
        if not ingested:
            return None
        return hduInfoList
