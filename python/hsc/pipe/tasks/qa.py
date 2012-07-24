from lsst.pex.config import Config, ConfigField
from lsst.pipe.base import Task
import hsc.pipe.tasks.measSeeingQa as hscSeeing

class QaConfig(Config):
    measureSeeing = ConfigField(dtype=hscSeeing.MeasureSeeingConfig, doc="Measure seeing")


class QaTask(Task):
    ConfigClass = QaConfig
    def __init__(self, *args, **kwargs):
        super(QaTask, self).__init(*args, **kwargs)
        self.makeSubtask("measureSeeing", hscSeeing.MeasureSeeingTask)

    def run(self, dataRef, exposure, sources):
        self.measureSeeing.run(dataRef, sources, exposure)

        metadata = exposure.getMetadata()

        # = flags
        # this part should be done by calculating merit functions somewhere else in a polite manner.
        metadata.set('FLAG_AUTO', 0)
        metadata.set('FLAG_USR', 0)
        metadata.set('FLAG_TAG', 1)

        # = info for data management
        # = rerun; assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'

        rerunName = self.getRerunName(dataRef)
        metadata.set('RERUN', rerunName)

    def getRerunName(self, sensorRef):
        return sensorRef.getButler().mapper.rerun

