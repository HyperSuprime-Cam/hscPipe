from lsst.meas.base.forcedPhotCcd import ForcedPhotCcdTask
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner
from lsst.pipe.base import ButlerInitializedTaskRunner

__all__ = ("ForcedPhotCcdDriverTask",)


class ForcedPhotCcdRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
    pass


class ForcedPhotCcdDriverTask(BatchParallelTask, ForcedPhotCcdTask):
    RunnerClass = ForcedPhotCcdRunner

