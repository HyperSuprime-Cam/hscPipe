from lsst.pex.config import Config, ConfigurableField
from lsst.pipe.base import Task
from . import measSeeingQa
from . import sizeMagnitudeMitakaStarSelector as starSel

class QaConfig(Config):
    ##seeing = ConfigurableField(target=measSeeingQa.MeasureSeeingTask, doc="Measure seeing")
    seeing = ConfigurableField(target=measSeeingQa.MeasureSeeingMitakaTask, doc="Measure seeing ver. Mitaka")

class QaTask(Task):
    ConfigClass = QaConfig

    def __init__(self, *args, **kwargs):
        super(QaTask, self).__init__(*args, **kwargs)
        self.makeSubtask("seeing", metadata=self.metadata)

        ## initialization of QA values
        # common
        self.metadata.set('RERUN', 'NOT_SET')
        self.metadata.set('FLAG_AUTO', -1)
        self.metadata.set('FLAG_USR', -1)
        self.metadata.set('FLAG_TAG', -1)
        # seeing/starsel
        self.metadata.set('SEEING_MODE', -1.0)
        self.metadata.set('ELL_MED', -1.0)
        self.metadata.set('ELL_PA_MED', -1.0)
        self.metadata.set("magLim", 99.0)
        self.metadata.set("fwhmRough", -1.0)
        self.metadata.set("fwhmRobust", -1.0)
        self.metadata.set("ellRobust", -1.0)
        self.metadata.set("ellPaRobust", -1.0)
        self.metadata.set("medianFwhmPsfSeq", -1.0)
        self.metadata.set("sigmaFwhmPsfSeq", -1.0)
        self.metadata.set("minFwhmPsfSeq", -1.0)
        self.metadata.set("maxFwhmPsfSeq", -1.0)
        self.metadata.set("numFwhmPsfLikeRobust", -1)

    def run(self, dataRef, exposure, sources):
        # initial set of exposure metadata
        #metadata = exposure.getMetadata()
        #metadata.combine(self.metadata)

        self.seeing.run(dataRef, sources, exposure)

        # = flags
        # this part should be done by calculating merit functions somewhere else in a polite manner.
        self.metadata.set('FLAG_AUTO', 0)
        self.metadata.set('FLAG_USR', 0)
        self.metadata.set('FLAG_TAG', 1)

        # = info for data management
        # = rerun; assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'
        try:
            rerunName = self.getRerunName(dataRef)
            self.metadata.set('RERUN', rerunName)
        except:
            self.log.warn("Could not determine rerun from output directory")

        # merging local Qa metadata into exposure metadata
        metadata = exposure.getMetadata()
        for key in self.metadata.names():
            metadata.set(key, self.metadata.get(key))
            #print '*** qa local metadata: %s = %s is set to exposure' % (key, str(self.metadata.get(key)))

    def getRerunName(self, sensorRef):
        # rerun: assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'
        corrPath = sensorRef.get('calexp_filename')[0]
        return corrPath[corrPath.find('rerun'):].split('/')[1]
