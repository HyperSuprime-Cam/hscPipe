import argparse
from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.tasks.forcedPhotCcd import ForcedPhotCcdTask
from lsst.pipe.tasks.dataIds import PerTractCcdDataIdContainer
from hsc.pipe.base.parallel import BatchPoolTask
from hsc.pipe.base.pool import Pool, abortOnError, NODE

class ForcedCcdConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    forcedPhotCcd = ConfigurableField(target=ForcedPhotCcdTask, doc="Forced measurement on ccd images")

    def validate(self):
        if self.forcedPhotCcd.references.coaddName != self.coaddName:
            raise RuntimeError("forcedPhotCcd.references.coaddName and coaddName don't match")

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

class ForcedCcdTask(BatchPoolTask):
    RunnerClass = ButlerInitializedTaskRunner
    ConfigClass = ForcedCcdConfig
    _DefaultName = "forcedCcd"

    def __init__(self, *args, **kwargs):
        try:
            butler = kwargs.pop("butler")
        except Exception, e:
            butler = None
        super(ForcedCcdTask, self).__init__(*args, **kwargs)
        if butler:
            self.makeSubtask("forcedPhotCcd", butler=butler)

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "calexp", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, dataRefList):
        """Run forced measurement on invididual ccd images

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param dataRefList: Data references to run measurement
        """
        pool = Pool("forcedCcd")
        ccdData = [dataRef for dataRef in dataRefList]
        pool.map(self.forced, ccdData)

    def forced(self, cache, dataRef):
        self.log.info("%s: Start forcedPhotCcd %s" % (NODE, dataRef.dataId))
        self.makeSubtask("forcedPhotCcd", butler=dataRef.getButler())
        self.forcedPhotCcd.run(dataRef)
        self.log.info("%s: Finished forcedPhotCcd %s" % (NODE, dataRef.dataId))

    def writeMetadata(self, dataRef):
        pass
