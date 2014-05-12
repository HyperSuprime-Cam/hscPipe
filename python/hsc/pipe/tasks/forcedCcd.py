from lsst.pipe.base import ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.tasks.forcedPhotCcd import ForcedPhotCcdTask, CcdForcedSrcDataIdContainer
from hsc.pipe.tasks.stack import TractDataIdContainer
from hsc.pipe.base.pbs import PbsPoolTask
from hsc.pipe.base.pool import Pool, abortOnError, NODE

class ForcedCcdDataIdContainer(CcdForcedSrcDataIdContainer):
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList
        """
        refList = dict()
        for dataId in self.idList:
            if "tract" not in dataId:
                raise argparse.ArgumentError(None, "--id must include tract")
            tract = dataId.pop("tract")
            if not tract in refList.keys():
                refList[tract] = list()
            # making a DataRef for src fills out any missing keys and allows us to iterate
            for srcDataRef in namespace.butler.subset("src", dataId=dataId):
                forcedDataId = srcDataRef.dataId.copy()
                forcedDataId['tract'] = tract
                dataRef = namespace.butler.dataRef(
                    datasetType = "forced_src",
                    dataId = forcedDataId,
                    )
                if not dataRef in refList[tract]:
                    refList[tract].append(dataRef)

        for tract in refList.keys():
            self.refList += [refList[tract]]

class ForcedCcdConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    forcedPhotCcd = ConfigurableField(target=ForcedPhotCcdTask, doc="Forced measurement on ccd images")

    def validate(self):
        if self.forcedPhotCcd.references.coaddName != self.coaddName:
            raise RuntimeError("forcedPhotCcd.references.coaddName and coaddName don't match")

class ForcedCcdTask(PbsPoolTask):
    ConfigClass = ForcedCcdConfig
    _DefaultName = "forcedCcd"

    def __init__(self, *args, **kwargs):
        super(ForcedCcdTask, self).__init__(*args, **kwargs)

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "forced_src", help="data ID, with raw CCD keys + tract",
                               ContainerClass=ForcedCcdDataIdContainer)
        return parser

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, patchRefList):
        print patchRefList
        """Run forced measurement on invididual ccd images

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param tractRef: Data reference for tract
        """
        pool = Pool("forcedCcd")
        coaddData = [patchRef for patchRef in patchRefList]
        pool.map(self.forced, coaddData)

    def forced(self, cache, patchRef):
        self.log.info("%s: Start forcedPhotCcd %s" % (NODE, patchRef.dataId))
        self.makeSubtask("forcedPhotCcd", butler=patchRef.getButler())
        self.forcedPhotCcd.run(patchRef)
        self.log.info("%s: Finished forcedPhotCcd %s" % (NODE, patchRef.dataId))

    def writeMetadata(self, dataRef):
        pass
