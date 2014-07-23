from lsst.pipe.base import ArgumentParser, TaskRunner
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

class ButlerInitializedTaskRunner(TaskRunner):
    """A TaskRunner for CmdLineTasks that require a 'butler' keyword argument to be passed to
    their constructor.
    """
    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor

        parsedCmd or args must be specified
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRef, kwargs = args
            butler = dataRef[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)

class ForcedCoaddTask(PbsPoolTask):
    RunnerClass = ButlerInitializedTaskRunner
    ConfigClass = ForcedCoaddConfig
    _DefaultName = "forcedCoadd"

    def __init__(self, *args, **kwargs):
        try:
            butler = kwargs.pop("butler")
        except Exception, e:
            butler = None
        super(ForcedCoaddTask, self).__init__(*args, **kwargs)
        if butler:
            self.makeSubtask("forcedPhotCoadd", butler=butler)

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
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

        @param patchRefList: Data references to run measurement
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
