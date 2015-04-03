from argparse import ArgumentError
from lsst.pex.config import Config, Field, ConfigurableField, FieldValidationError
from lsst.pipe.base.argumentParser import ArgumentParser
from lsst.pipe.tasks.multiBand import (DetectCoaddSourcesTask, MergeDetectionsTask,
                                       MeasureMergedCoaddSourcesTask, MergeMeasurementsTask,)
from lsst.pipe.tasks.forcedPhotCoadd import ForcedPhotCoaddTask
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, CoaddTaskRunner
from lsst.pipe.tasks.references import CoaddSrcReferencesTask, MultiBandReferencesTask
from hsc.pipe.tasks.stack import TractDataIdContainer
from hsc.pipe.base.parallel import BatchPoolTask
from hsc.pipe.base.pool import Pool, abortOnError
from hsc.pipe.base.butler import getDataRef

import lsst.afw.table as afwTable

class MultiBandDataIdContainer(CoaddDataIdContainer):
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList

        It's difficult to make a data reference that merely points to an entire
        tract: there is no data product solely at the tract level.  Instead, we
        generate a list of data references for patches within the tract.
        """
        datasetType = namespace.config.coaddName + "Coadd"
        getPatchRefList = lambda tract: [namespace.butler.dataRef(datasetType=datasetType,
                                                                  tract=tract.getId(),
                                                                  filter=dataId["filter"],
                                                                  patch="%d,%d" % patch.getIndex()) for
                                         patch in tract]

        tractRefs = {} # Data references for each tract
        for dataId in self.idList:
            # There's no registry of coadds by filter, so we need to be given the filter
            if "filter" not in dataId:
                raise ArgumentError(None, "--id must include 'filter'")

            skymap = self.getSkymap(namespace, datasetType)

            if "tract" in dataId:
                tractId = dataId["tract"]
                if tractId not in tractRefs:
                    tractRefs[tractId] = []
                if "patch" in dataId:
                    tractRefs[tractId].append(namespace.butler.dataRef(datasetType=datasetType,
                                                                       tract=tractId,
                                                                       filter=dataId['filter'],
                                                                       patch=dataId['patch']))
                else:
                    tractRefs[tractId] += getPatchRefList(skymap[tractId])
            else:
                tractRefs = dict((tract.getId(), tractRefs.get(tract.getId(), []) + getPatchRefList(tract))
                                 for tract in skymap)

        self.refList = tractRefs.values()


class MultiBandTaskRunner(CoaddTaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        kwargs["butler"] = parsedCmd.butler
        return [(parsedCmd.id.refList, kwargs),]


class MultiBandConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    detectCoaddSources = ConfigurableField(target=DetectCoaddSourcesTask, doc="Detect sources on coadd")
    mergeCoaddDetections = ConfigurableField(target=MergeDetectionsTask, doc="Merge detections")
    measureCoaddSources = ConfigurableField(target=MeasureMergedCoaddSourcesTask,
                                            doc="Measure merged detections")
    mergeCoaddMeasurements = ConfigurableField(target=MergeMeasurementsTask, doc="Merge measurements")
    forcedPhotCoadd = ConfigurableField(target=ForcedPhotCoaddTask,
                                        doc="Forced measurement on coadded images")

    def setDefaults(self):
        Config.setDefaults(self)
        self.forcedPhotCoadd.references.retarget(MultiBandReferencesTask)

    def validate(self):
        for subtask in ("detectCoaddSources", "mergeCoaddDetections", "measureCoaddSources",
                        "mergeCoaddMeasurements", "forcedPhotCoadd"):
            coaddName = getattr(self, subtask).coaddName
            if coaddName != self.coaddName:
                raise RuntimeError("%s.coaddName (%s) doesn't match root coaddName (%s)" %
                                   (subtask, coaddName, self.coaddName))

class MultiBandTask(BatchPoolTask):
    """Multi-node driver for multiband processing"""
    ConfigClass = MultiBandConfig
    _DefaultName = "multiband"

    def __init__(self, butler=None, **kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        self.makeSubtask("detectCoaddSources", butler=butler) # Schema may have been written in stack.py
        self.makeSubtask("mergeCoaddDetections", schema=afwTable.Schema(self.detectCoaddSources.schema))
        self.makeSubtask("measureCoaddSources", butler=butler,
                         schema=afwTable.Schema(self.mergeCoaddDetections.schema),
                         peakSchema=afwTable.Schema(self.mergeCoaddDetections.merged.getPeakSchema()))
        self.makeSubtask("mergeCoaddMeasurements", butler=butler,
                         schema=afwTable.Schema(self.measureCoaddSources.schema))
        self.makeSubtask("forcedPhotCoadd", butler=butler,
                         schema=afwTable.Schema(self.mergeCoaddMeasurements.schema))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="multiband", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = 0
        for refList in parsedCmd.id.refList:
            numTargets += len(refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, patchRefList):
        """Run multiband processing on coadds

        All nodes execute this method, though the master and slaves
        take different routes through it.

        No real MPI communication takes place: all I/O goes through the disk.
        We want the intermediate stages on disk, and the component Tasks are
        implemented around this, so we just follow suit.

        patchRefList:  Data references to run measurement
        """
        for patchRef in patchRefList:
            if patchRef:
                butler = patchRef.getButler()
                break
        else:
            raise RuntimeError("No valid patches")
        pool = Pool("all")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        patchRefList = [patchRef for patchRef in patchRefList if
                           patchRef.datasetExists(self.config.coaddName + "Coadd")]
        dataIdList = [patchRef.dataId for patchRef in patchRefList]

        # Group by patch
        patches = {}
        tract = None
        for patchRef in patchRefList:
            dataId = patchRef.dataId
            if tract is None:
                tract = dataId["tract"]
            else:
                assert tract == dataId["tract"]

            patch = dataId["patch"]
            if patch not in patches:
                patches[patch] = []
            patches[patch].append(dataId)

        pool.map(self.runDetect, dataIdList)
        pool.map(self.runMergeDetections, patches.values())
        pool.map(self.runMeasureMerged, dataIdList)
        pool.map(self.runMergeMeasurements, patches.values())
        pool.map(self.runForcedPhot, dataIdList)

    def runDetect(self, cache, dataId):
        """Run detection on a patch for a single filter

        Only slave nodes execute this method.

        cache: Pool cache, containing butler
        dataId: Data identifier for patch
        """
        with self.logOperation("detection on %s" % (dataId,)):
            dataRef = getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd")
            if dataRef.datasetExists(self.config.coaddName + "Coadd_det"):
                return
            self.detectCoaddSources.run(dataRef)

    def runMergeDetections(self, cache, dataIdList):
        """Run detection merging on a patch

        Only slave nodes execute this method.

        cache: Pool cache, containing butler
        dataIdList: List of data identifiers for the patch in different filters
        """
        with self.logOperation("merge detections from %s" % (dataIdList,)):
            dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd") for
                           dataId in dataIdList]
            if dataRefList[0].datasetExists(self.config.coaddName + "Coadd_mergeDet"):
                return
            self.mergeCoaddDetections.run(dataRefList)

    def runMeasureMerged(self, cache, dataId):
        """Run measurement on a patch for a single filter

        Only slave nodes execute this method.

        cache: Pool cache, with butler
        dataId: Data identifier for patch
        """
        with self.logOperation("measurement on %s" % (dataId,)):
            dataRef = getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd")
            if dataRef.datasetExists(self.config.coaddName + "Coadd_meas"):
                return
            self.measureCoaddSources.run(dataRef)

    def runMergeMeasurements(self, cache, dataIdList):
        """Run measurement merging on a patch

        Only slave nodes execute this method.

        cache: Pool cache, containing butler
        dataIdList: List of data identifiers for the patch in different filters
        """
        with self.logOperation("merge measurements from %s" % (dataIdList,)):
            dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd") for
                           dataId in dataIdList]
            if dataRefList[0].datasetExists(self.config.coaddName + "Coadd_ref"):
                return
            self.mergeCoaddMeasurements.run(dataRefList)

    def runForcedPhot(self, cache, dataId):
        """Run forced photometry on a patch for a single filter

        Only slave nodes execute this method.

        cache: Pool cache, with butler
        dataId: Data identifier for patch
        """
        with self.logOperation("forced photometry on %s" % (dataId,)):
            dataRef = getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd")
            if dataRef.datasetExists(self.config.coaddName + "Coadd_forced_src"):
                return
            self.forcedPhotCoadd.run(dataRef)

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass
