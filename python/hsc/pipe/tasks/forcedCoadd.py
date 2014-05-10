from lsst.pipe.base import ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.tasks.forcedPhotCoadd import ForcedPhotCoaddTask
from hsc.pipe.tasks.stack import TractDataIdContainer
from hsc.pipe.base.pbs import PbsPoolTask
from hsc.pipe.base.pool import Pool, abortOnError, NODE

class ForcedCoaddConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    forcedPhotCoadd = ConfigurableField(target=ForcedPhotCoaddTask, doc="Forced measurement on coadded images")

    def validate(self):
        if self.forcedPhotCoadd.coaddName != self.coaddName:
            raise RuntimeError("forcedPhotCoadd.coaddName and coaddName don't match")

class ForcedCoaddTask(PbsPoolTask):
    ConfigClass = ForcedCoaddConfig
    _DefaultName = "forcedCoadd"

    def __init__(self, *args, **kwargs):
        super(ForcedCoaddTask, self).__init__(*args, **kwargs)

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
        """
        Selection references are not cheap (reads Wcs), so are generated
        only if we're not doing a PBS submission.
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, patchRefList):
        """Run forced measurement on coadding images

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param tractRef: Data reference for tract
        """
        pool = Pool("forcedCoadd")
        coaddData = [patchRef for patchRef in patchRefList]
        pool.map(self.forced, coaddData)

    def forced(self, cache, patchRef):
        self.log.info("%s: Start forcedPhotCoadd %s" % (NODE, patchRef.dataId))
        self.makeSubtask("forcedPhotCoadd", butler=patchRef.getButler())
        self.forcedPhotCoadd.run(patchRef)
        self.log.info("%s: Finished forcedPhotCoadd %s" % (NODE, patchRef.dataId))

    def writeMetadata(self, dataRef):
        pass
