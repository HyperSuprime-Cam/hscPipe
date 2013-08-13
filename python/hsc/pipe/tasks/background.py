import argparse
from lsst.pex.config import Config, ConfigurableField, Field
from lsst.pipe.base import Task, CmdLineTask, Struct, ArgumentParser
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner, CoaddDataIdContainer, SelectDataIdContainer
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
import lsst.afw.geom as afwGeom

class AssignTask(Task):
    ConfigClass = Config

    def run(self, tractInfo, dataRefList):
        # XXX placeholder
        return None


class TractDataIdContainer(CoaddDataIdContainer):
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList"""
        datasetType = namespace.config.coaddName + "Coadd_bgRef"
        validKeys = namespace.butler.getKeys(datasetType=datasetType, level=self.level)

        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract",):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(None, "--id must include " + key)

            # tract is required; iterate over it if not provided
            if not "tract" in dataId:
                addList = [dict(tract=tract.getId(), **dataId) for
                           tract in self.getSkymap(namespace, datasetType)]
            else:
                addList = [dataId]

            self.refList += [namespace.butler.dataRef(datasetType=datasetType, dataId=addId)
                             for addId in addList]

class BackgroundReferenceIoTask(Task):
    """Provides abstraction of I/O for background references

    This implementation uses pickle files; one could imagine using
    a database instead.
    """
    ConfigClass = Config

    def __init__(self, *args, **kwargs):
        super(BackgroundReferenceIoTask, self).__init__(*args, **kwargs)
        import cPickle as pickle
        self._pickle = pickle

    @classmethod
    def _bgRefName(cls, coaddName):
        return coaddName + "Coadd_bgRef"

    @classmethod
    def _filename(cls, tractRef, coaddName):
        return tractRef.get(cls._bgRefName(coaddName) + "_filename", immediate=True)[0]

    def exists(self, tractRef, coaddName):
        return tractRef.datasetExists(self._bgRefName(coaddName))

    def write(self, tractRef, coaddName, assignments):
        f = open(self._filename(tractRef, coaddName), "w")
        self._pickle.dump(assignments, f)

    def read(self, tractRef, coaddName):
        f = open(self._filename(tractRef, coaddName))
        return self._pickle.load(f)


class BackgroundReferenceConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd of interest")
    io = ConfigurableField(target=BackgroundReferenceIoTask, doc="I/O for background references")
    select = ConfigurableField(target=WcsSelectImagesTask, doc="Task to select input images")
    assign = ConfigurableField(target=AssignTask, doc="Task to assign background refernces")
    clobber = Field(dtype=bool, default=False, doc="Clobber existing assignments?")

class BackgroundReferenceTask(CmdLineTask):
    _DefaultName = "bgRef"
    ConfigClass = BackgroundReferenceConfig
    RunnerClass = CoaddTaskRunner

    def __init__(self, config, **kwargs):
        super(BackgroundReferenceTask, self).__init__(config=config, **kwargs)
        self.makeSubtask("io")
        self.makeSubtask("select")
        self.makeSubtask("assign")
        self.bgRefName = self.config.coaddName + "Coadd_bgRef"

    @classmethod
    def _makeArgumentParser(cls):
        """Create an argument parser
        """
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345",
                               ContainerClass=TractDataIdContainer)
        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectDataIdContainer)
        return parser

    def run(self, tractRef, selectDataList=[]):
        if self.io.exists(tractRef, self.config.coaddName):
            if not self.config.clobber:
                self.log.warn("Refusing to clobber existing background reference for %s" % (tractRef,))
                return
            self.log.warn("Clobbering existing background reference for %s" % (tractRef,))

        skyMap = tractRef.get(self.config.coaddName + "Coadd_skyMap")
        tract = skyMap[tractRef.dataId['tract']]
        dataRefList = self.selectExposures(tractRef, selectDataList, tractInfo=tract)
        assignments = self.assign.run(tract, dataRefList)
        self.io.write(tractRef, self.config.coaddName, assignments)

    def selectExposures(self, tractRef, selectDataList=[], tractInfo=None):
        """Select exposures to include

        @param tractRef: data reference for sky map patch
        @param selectDataList: list of SelectStruct
        @param tractInfo: tract object from skymap
        @return list of calexp data references
        """
        if tractInfo is None:
            skyMap = tractRef.get(coaddName + "Coadd_skyMap")
            tractInfo = skyMap[tractRef.dataId['tract']]
        wcs = tractInfo.getWcs()
        cornerPosList = afwGeom.Box2D(tractInfo.getBBox()).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        return self.select.runDataRef(tractRef, coordList, selectDataList=selectDataList).dataRefList

    def writeMetadata(self, *args, **kwargs):
        pass
    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchemas(self, *args, **kwargs):
        pass
