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

class ForcedCcdTaskRunner(TaskRunner):
    """Provide a list of data references as targets instead of just one at a time

    This allows ForcedCcdTask.run to iterate over them using the Pool.
    """
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        return [(parsedCmd.id.refList, kwargs)]

    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor

        parsedCmd or args must be specified

        This is like the ButlerInitializedTaskRunner except the 'args' contains
        a list of data references rather than a single data reference.
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)

class ForcedCcdTask(BatchPoolTask):
    RunnerClass = ForcedCcdTaskRunner
    ConfigClass = ForcedCcdConfig
    _DefaultName = "forcedCcd"

    def __init__(self, butler, *args, **kwargs):
        super(ForcedCcdTask, self).__init__(*args, **kwargs)
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
        numTargets = len(parsedCmd.id.refList)
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
        self.forcedPhotCcd.run(dataRef)
        self.log.info("%s: Finished forcedPhotCcd %s" % (NODE, dataRef.dataId))

    def writeMetadata(self, dataRef):
        pass
