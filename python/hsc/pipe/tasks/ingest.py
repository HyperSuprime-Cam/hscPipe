from glob import glob
from lsst.pex.config import Config, ConfigurableField
from lsst.obs.subaru.ingest import HscIngestTask, HscIngestArgumentParser
from hsc.pipe.base.pool import Pool, startPool, abortOnError
from hsc.pipe.base.pbs import PbsCmdLineTask

class PoolIngestTask(PbsCmdLineTask, HscIngestTask):
    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        return float(time)*len(parsedCmd.files)/numNodes

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        return HscIngestArgumentParser("ingest")

    @classmethod
    def parseAndRun(cls, *args, **kwargs):
        """Run with a MPI process pool"""
        pool = startPool()
        config = cls.ConfigClass()
        parser = cls.ArgumentParser("ingest")
        args = parser.parse_args(config)
        task = cls(config=args.config)
        task.run(args)
        pool.exit()

    @abortOnError
    def run(self, args):
        pool = Pool(None)

        # Read and move files in parallel
        filenameList = sum([glob(filename) for filename in args.files], [])
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
