import os
from argparse import ArgumentError
from lsst.pex.config import Config, Field, ConfigurableField, FieldValidationError
from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.pipe.tasks.multiBand import (MergeDetectionsTask,
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
        datasetType = namespace.config.coaddName + "Coadd_calexp"
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
    mergeCoaddDetections = ConfigurableField(target=MergeDetectionsTask, doc="Merge detections")
    measureCoaddSources = ConfigurableField(target=MeasureMergedCoaddSourcesTask,
                                            doc="Measure merged detections")
    mergeCoaddMeasurements = ConfigurableField(target=MergeMeasurementsTask, doc="Merge measurements")
    forcedPhotCoadd = ConfigurableField(target=ForcedPhotCoaddTask,
                                        doc="Forced measurement on coadded images")
    clobberDetections = Field(dtype=bool, default=False, doc="Clobber existing detections?")
    clobberMergedDetections = Field(dtype=bool, default=False, doc="Clobber existing merged detections?")
    clobberMeasurements = Field(dtype=bool, default=False, doc="Clobber existing measurements?")
    clobberMergedMeasurements = Field(dtype=bool, default=False, doc="Clobber existing merged measurements?")
    clobberForcedPhotometry = Field(dtype=bool, default=False, doc="Clobber existing forced photometry?")
    reprocessing = Field(
        dtype=bool, default=False,
        doc=("Are we reprocessing?\n\n"
             "This exists as a workaround for large deblender footprints causing large memory use and/or very "
             "slow processing.  We refuse to deblend those footprints when running on a cluster and return to "
             "reprocess on a machine with larger memory or more time if we consider those footprints "
             "important to recover."),
        )

    def setDefaults(self):
        Config.setDefaults(self)
        self.forcedPhotCoadd.references.retarget(MultiBandReferencesTask)

    def validate(self):
        for subtask in ("mergeCoaddDetections", "measureCoaddSources",
                        "mergeCoaddMeasurements", "forcedPhotCoadd"):
            coaddName = getattr(self, subtask).coaddName
            if coaddName != self.coaddName:
                raise RuntimeError("%s.coaddName (%s) doesn't match root coaddName (%s)" %
                                   (subtask, coaddName, self.coaddName))

class MultiBandTaskRunner(TaskRunner):
    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor

        parsedCmd or args must be specified
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)

def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class MultiBandTask(BatchPoolTask):
    """Multi-node driver for multiband processing"""
    ConfigClass = MultiBandConfig
    _DefaultName = "multiband"
    RunnerClass = MultiBandTaskRunner

    def __init__(self, butler=None, schema=None, **kwargs):
        BatchPoolTask.__init__(self, **kwargs)
        if schema is None:
            assert butler is not None, "Butler not provided"
            schema = butler.get(self.config.coaddName + "Coadd_det_schema", immediate=True).schema
        self.butler = butler
        self.makeSubtask("mergeCoaddDetections", schema=schema)
        self.makeSubtask("measureCoaddSources", schema=afwTable.Schema(self.mergeCoaddDetections.schema),
                         peakSchema=afwTable.Schema(self.mergeCoaddDetections.merged.getPeakSchema()))
        self.makeSubtask("mergeCoaddMeasurements", schema=afwTable.Schema(self.measureCoaddSources.schema))
        self.makeSubtask("forcedPhotCoadd", schema=afwTable.Schema(self.mergeCoaddMeasurements.schema))

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   butler=self.butler))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="multiband", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd_calexp", help="data ID, e.g. --id tract=12345 patch=1,2",
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
                        patchRef.datasetExists(self.config.coaddName + "Coadd_calexp")]
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

        pool.map(self.runMergeDetections, patches.values())

        # Measure merged detections, and test for reprocessing
        #
        # The reprocessing allows us to have multiple attempts at deblending large footprints. Large
        # footprints can suck up a lot of memory in the deblender, which means that when we process on a
        # cluster, we want to refuse to deblend them (they're flagged "deblend.parent-too-big"). But since
        # they may have astronomically interesting data, we want the ability to go back and reprocess them
        # with a more permissive configuration when we have more memory or processing time.
        #
        # self.runMeasureMerged will return whether there are any footprints in that image that required
        # reprocessing.  We need to convert that list of booleans into a dict mapping the patchId (x,y) to
        # a boolean. That tells us whether the merge measurement and forced photometry need to be re-run on
        # a particular patch.
        #
        # This determination of which patches need to be reprocessed exists only in memory (the measurements
        # have been written, clobbering the old ones), so if there was an exception we would lose this
        # information, leaving things in an inconsistent state (measurements new, but merged measurements and
        # forced photometry old). To attempt to preserve this status, we touch a file (dataset named
        # "deepCoadd_multibandReprocessing") --- if this file exists, we need to re-run the merge and
        # forced photometry.
        #
        # This is, hopefully, a temporary workaround until we can improve the deblender.
        try:
            reprocessed = pool.map(self.runMeasureMerged, dataIdList)
        finally:
            if self.config.reprocessing:
                patchReprocessing = {}
                for dataId, reprocess in zip(dataIdList, reprocessed):
                    patchId = dataId["patch"]
                    patchReprocessing[patchId] = patchReprocessing.get(patchId, False) or reprocess
                # Persist the determination, to make error recover easier
                reprocessDataset = self.config.coaddName + "Coadd_multibandReprocessing"
                for patchId in patchReprocessing:
                    if not patchReprocessing[patchId]:
                        continue
                    dataId = dict(tract=tract, patch=patchId)
                    if patchReprocessing[patchId]:
                        filename = butler.get(reprocessDataset + "_filename", dataId)[0]
                        open(filename, 'a').close() # Touch file
                    elif butler.datasetExists(reprocessDataset, dataId):
                        # We must have failed at some point while reprocessing and we're starting over
                        patchReprocessing[patchId] = True

        # Only process patches that have been identified as needing it
        pool.map(self.runMergeMeasurements, [idList for patchId, idList in patches.iteritems() if
                                             not self.config.reprocessing or patchReprocessing[patchId]])
        pool.map(self.runForcedPhot, [dataId for dataId in dataIdList if not self.config.reprocessing or
                                      patchReprocessing[dataId["patch"]]])

        # Remove persisted reprocessing determination
        if self.config.reprocessing:
            for patchId in patchReprocessing:
                if not patchReprocessing[patchId]:
                    continue
                dataId = dict(tract=tract, patch=patchId)
                filename = butler.get(reprocessDataset + "_filename", dataId)[0]
                os.unlink(filename)

    def runMergeDetections(self, cache, dataIdList):
        """Run detection merging on a patch

        Only slave nodes execute this method.

        cache: Pool cache, containing butler
        dataIdList: List of data identifiers for the patch in different filters
        """
        with self.logOperation("merge detections from %s" % (dataIdList,)):
            dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                           dataId in dataIdList]
            if (not self.config.clobberMergedDetections and
                dataRefList[0].datasetExists(self.config.coaddName + "Coadd_mergeDet")):
                return
            self.mergeCoaddDetections.run(dataRefList)

    def runMeasureMerged(self, cache, dataId):
        """Run measurement on a patch for a single filter

        Only slave nodes execute this method.

        cache: Pool cache, with butler
        dataId: Data identifier for patch

        Returns whether the patch requires reprocessing.
        """
        with self.logOperation("measurement on %s" % (dataId,)):
            dataRef = getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp")
            reprocessing = False # Does this patch require reprocessing?
            if (not self.config.clobberMeasurements and
                dataRef.datasetExists(self.config.coaddName + "Coadd_meas")):
                if not self.config.reprocessing:
                    return False

                catalog = dataRef.get(self.config.coaddName + "Coadd_meas")
                bigFlag = catalog["deblend.parent-too-big"]
                numOldBig = bigFlag.sum()
                if numOldBig == 0:
                    self.log.info("No large footprints in %s" % (dataRef.dataId,))
                    return False
                numNewBig = sum((self.measureCoaddSources.deblend.isLargeFootprint(src.getFootprint()) for
                                 src in catalog[bigFlag]))
                if numNewBig == numOldBig:
                    self.log.info("All %d formerly large footprints continue to be large in %s" %
                                  (numOldBig, dataRef.dataId,))
                    return False
                self.log.info("Found %d large footprints to be reprocessed in %s" %
                              (numOldBig - numNewBig, dataRef.dataId))
                reprocessing = True

            self.measureCoaddSources.run(dataRef)
            return reprocessing

    def runMergeMeasurements(self, cache, dataIdList):
        """Run measurement merging on a patch

        Only slave nodes execute this method.

        cache: Pool cache, containing butler
        dataIdList: List of data identifiers for the patch in different filters
        """
        with self.logOperation("merge measurements from %s" % (dataIdList,)):
            dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                           dataId in dataIdList]
            if (not self.config.clobberMergedMeasurements and
                not self.config.reprocessing and
                dataRefList[0].datasetExists(self.config.coaddName + "Coadd_ref")):
                return
            self.mergeCoaddMeasurements.run(dataRefList)

    def runForcedPhot(self, cache, dataId):
        """Run forced photometry on a patch for a single filter

        Only slave nodes execute this method.

        cache: Pool cache, with butler
        dataId: Data identifier for patch
        """
        with self.logOperation("forced photometry on %s" % (dataId,)):
            dataRef = getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp")
            if (not self.config.clobberForcedPhotometry and
                not self.config.reprocessing and
                dataRef.datasetExists(self.config.coaddName + "Coadd_forced_src")):
                return
            self.forcedPhotCoadd.run(dataRef)

    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass
