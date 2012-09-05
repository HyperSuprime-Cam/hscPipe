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
        self.metadata.set('REGISTID', -9999)
        self.metadata.set('VISIT', -9999)
        self.metadata.set('CCD_REGISTRY', -9999)
        self.metadata.set('ANAPATH', 'NOT_SET')
        self.metadata.set('RERUN', 'NOT_SET')
        self.metadata.set('FLAG_AUTO', -9999)
        self.metadata.set('FLAG_USR', -9999)
        self.metadata.set('FLAG_TAG', -9999)
        # seeing/starsel
        self.metadata.set('SEEING_MODE', -9999.0)
        self.metadata.set('ELL_MED', -9999.0)
        self.metadata.set('ELL_PA_MED', -9999.0)
        self.metadata.set("magLim", 99.0)
        self.metadata.set("fwhmRough", -9999.0)
        self.metadata.set("fwhmRobust", -9999.0)
        self.metadata.set("ellRobust", -9999.0)
        self.metadata.set("ellPaRobust", -9999.0)
        self.metadata.set("medianFwhmPsfSeq", -9999.0)
        self.metadata.set("sigmaFwhmPsfSeq", -9999.0)
        self.metadata.set("minFwhmPsfSeq", -9999.0)
        self.metadata.set("maxFwhmPsfSeq", -9999.0)
        self.metadata.set("numFwhmPsfLikeRobust", -9999)

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

        # registry information
        try:
            anaPath = dataRef.getButler().mapper.root
            dataId = dataRef.dataId
            if dataId.has_key('registryId'):
                registryId = dataId['registryId']
            else:
                registryId = -9999
            visit = dataId['visit']
            ccd = dataId['ccd']
            self.metadata.set('REGISTID', registryId)
            self.metadata.set('VISIT', visit)
            self.metadata.set('CCD_REGISTRY', ccd)
            self.metadata.set('ANAPATH', anaPath)
            
            print '*** registry info dataRef.dataId: %s (visit, ccd, registryId): %s' % (str(dataId), str((visit, ccd, registryId)))

        except Exception, e:
            self.log.warn("Could not get registry info (visit, ccd): %s" % str(e))

        # rerun; assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'
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
