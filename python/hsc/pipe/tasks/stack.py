import re
import argparse

import lsst.afw.geom as afwGeom
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pex.exceptions import LsstCppException, InvalidParameterException
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, SelectDataIdContainer, CoaddTaskRunner
from lsst.pipe.tasks.selectImages import BaseSelectImagesTask, WcsSelectImagesTask, BaseExposureInfo
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.tasks.assembleCoadd import AssembleCoaddTask
from lsst.pipe.tasks.processCoadd import ProcessCoaddTask
from lsst.pipe.base import Struct, DataIdContainer
from hsc.pipe.base.pbs import PbsCmdLineTask
from hsc.pipe.base.mpi import (MpiTask, MpiMultiplexTaskRunner, MpiArgumentParser, thisNode, abortOnError)
from hsc.pipe.tasks.background import MpiBackgroundReferenceTask

###for cls in (MakeCoaddTempExpTask, AssembleCoaddTask):
###    cls = wrapTask(cls, globals())
###
###    # Replace the _makeArgumentParser classmethod with one that's MPI-friendly.
###    # The result is that the selection data is only defined on the root node.
###    # This is important because we don't want all nodes reading the WCSes of all the inputs at once.
###    def _makeArgumentParser(cls):
###        parser = MpiArgumentParser(name=cls._DefaultName)
###        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
###                               ContainerClass=CoaddDataIdContainer, rootOnly=False)
###        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
###                               ContainerClass=SelectDataIdContainer, rootOnly=True)
###        return parser
###    cls._makeArgumentParser = types.MethodType(_makeArgumentParser, cls, cls.__metaclass__)
###
###    # Use a SelectImagesTask that will transfer the selection data over MPI.
###    class WrapConfig(cls.ConfigClass):
###        def setDefaults(self):
###            self.selectImages.retarget(MpiWcsSelectImagesTask)
###
###    cls.ConfigClass = WrapConfig


class NullSelectImagesTask(BaseSelectImagesTask):
    """Select images by taking everything we're given without further examination

    This is useful if the examination (e.g., Wcs checking) has been performed
    previously, and we've been provided a good list.
    """
    def runDataRef(self, patchRef, coordList, makeDataRefList=True, selectDataList=[]):
        return Struct(
            dataRefList = [s.dataRef for s in selectDataList],
            exposureInfoList = [BaseExposureInfo(s.dataRef.dataId, None) for s in selectDataList],
            )


class StackConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    select = ConfigurableField(target=BaseSelectImagesTask, doc="Select images to process")
    makeCoaddTempExp = ConfigurableField(target=MakeCoaddTempExpTask, doc="Warp images to sky")
    backgroundReference = ConfigurableField(target=MpiBackgroundReferenceTask,
                                            doc="Build background reference")
    assembleCoadd = ConfigurableField(target=AssembleCoaddTask, doc="Assemble warps into coadd")
    processCoadd = ConfigurableField(target=ProcessCoaddTask, doc="Detection and measurement on coadd")
    doOverwriteCoadd = Field(dtype=bool, default=False, doc="Overwrite coadd?")
    doOverwriteOutput = Field(dtype=bool, default=False, doc="Overwrite processing outputs?")

    def setDefaults(self):
        self.select.retarget(WcsSelectImagesTask)
        self.makeCoaddTempExp.select.retarget(NullSelectImagesTask)
        self.backgroundReference.construct.doWarp = False
        self.backgroundReference.select.retarget(NullSelectImagesTask)
        self.assembleCoadd.select.retarget(NullSelectImagesTask)
        self.makeCoaddTempExp.doOverwrite = False

        # Not quite ready for background matching yet
        self.assembleCoadd.doMatchBackgrounds = False
        self.makeCoaddTempExp.bgSubtracted = True

    def validate(self):
        if self.makeCoaddTempExp.coaddName != self.coaddName:
            raise RuntimeError("makeCoaddTempExp.coaddName and coaddName don't match")
        if self.assembleCoadd.coaddName != self.coaddName:
            raise RuntimeError("assembleCoadd.coaddName and coaddName don't match")
        if self.backgroundReference.coaddName != self.coaddName:
            raise RuntimeError("backgroundReference.coaddName and coaddName don't match")

class TractDataIdContainer(CoaddDataIdContainer):
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList

        It's difficult to make a data reference that merely points to an entire
        tract: there is no data product solely at the tract level.  Instead, we
        generate a list of data references for patches within the tract.
        """
        datasetType = namespace.config.coaddName + "Coadd"
        validKeys = set(["tract", "filter"])

        getPatchRefList = lambda tract: [namespace.butler.dataRef(datasetType=datasetType, tract=tract.getId(),
                                                                  filter=dataId["filter"],
                                                                  patch="%d,%d" % patch.getIndex()) for
                                         patch in tract]

        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract"):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(None, "--id must include " + key)

            skymap = self.getSkymap(namespace, datasetType)

            # tract is required; iterate over it if not provided
            if not "tract" in dataId:
                addList = [getPatchRefList(tract) for tract in skymap]
            else:
                addList = [getPatchRefList(skymap[dataId["tract"]])]
            self.refList += addList

class StackTask(PbsCmdLineTask, MpiTask):
    ConfigClass = StackConfig
    _DefaultName = "stacker" # "stack" conflicts with hscMosaic's StackTask.
    RunnerClass = CoaddTaskRunner

    def __init__(self, *args, **kwargs):
        super(StackTask, self).__init__(*args, **kwargs)
        self.makeSubtask("select")
        self.makeSubtask("makeCoaddTempExp")
        self.makeSubtask("backgroundReference")
        self.makeSubtask("assembleCoadd")
        self.makeSubtask("processCoadd")

    @classmethod
    def _makeArgumentParser(cls, doPbs=False, **kwargs):
        """
        Patch references are cheap, so are defined on all nodes.
        Selection references are not cheap (reads Wcs), so are generated
        only on the root node (and only if we're not doing a PBS submission),
        and distributed via the MpiWcsSelectImagesTask.
        """
        parser = MpiArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer, rootOnly=False)
        # We don't want to be reading all the WCSes if we're only in the act of submitting to PBS
        SelectContainerClass = DataIdContainer if doPbs else SelectDataIdContainer
        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectContainerClass, rootOnly=True)
        return parser

    @classmethod
    def pbsWallTime(cls, time, parsedCmd, numNodes, numProcs):
        numTargets = len(parsedCmd.selectId.refList)
        return time*numTargets/float(numNodes*numProcs)

    @abortOnError
    def run(self, patchRefList, selectDataList=[]):
        """Run stacking on a tract

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param tractRef: Data reference for tract
        @param selectDataList: List of SelectStruct for inputs
        """
        import pbasf2
        if self.rank == self.root:
            warpData = [Struct(patchRef=patchRef, selectDataList=selectDataList) for patchRef in patchRefList]
        else:
            warpData = None
        selectedData = pbasf2.ScatterJob(self.comm, self.warp, warpData, root=self.root)
#        self.backgroundReference.run(patchRefList, selectDataList)

        if self.rank == self.root:
            refNamer = lambda patchRef: tuple(map(int, patchRef.dataId["patch"].split(",")))
            lookup = dict(zip(map(refNamer, patchRefList), selectedData))
            coaddData = [Struct(patchRef=patchRef, selectDataList=lookup[refNamer(patchRef)]) for
                         patchRef in patchRefList]
        else:
            coaddData = None
        pbasf2.ScatterJob(self.comm, self.coadd, coaddData, root=self.root)

    def warp(self, data):
        """Warp all images for a patch

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param selectDataList: List of SelectStruct for inputs
        @return selectDataList with non-overlapping elements removed
        """
        patchRef = data.patchRef
        selectDataList = data.selectDataList
        self.log.info("%s: Start warping %s" % (thisNode(), patchRef.dataId))
        selectDataList = self.selectExposures(patchRef, selectDataList)
        self.makeCoaddTempExp.run(patchRef, selectDataList)
        self.log.info("%s: Finished warping %s" % (thisNode(), patchRef.dataId))
        return selectDataList

    def coadd(self, data):
        """Construct coadd for a patch and measure

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param selectDataList: List of SelectStruct for inputs
        """
        patchRef = data.patchRef
        selectDataList = data.selectDataList
        self.log.info("%s: Start coadding %s" % (thisNode(), patchRef.dataId))
        coaddName = self.config.coaddName + "Coadd"
        coadd = None
        if self.config.doOverwriteCoadd or not patchRef.datasetExists(coaddName):
            coaddResults = self.assembleCoadd.run(patchRef, selectDataList)
            if coaddResults is not None:
                coadd = coaddResults.coaddExposure
        elif patchRef.datasetExists(coaddName):
            coadd = patchRef.get(coaddName, immediate=True)
        self.log.info("%s: Finished coadding %s" % (thisNode(), patchRef.dataId))
        if coadd is not None and (self.config.doOverwriteOutput or
                                  not patchRef.datasetExists(coaddName + "_src") or
                                  not patchRef.datasetExists(coaddName + "_calexp")):
            self.process(patchRef, coadd)

    def selectExposures(self, patchRef, selectDataList):
        """Select exposures to operate upon, via the SelectImagesTask

        This is very similar to CoaddBaseTask.selectExposures, except we return
        a list of SelectStruct (same as the input), so we can plug the results into
        future uses of SelectImagesTask.
        """
        key = lambda dataRef: tuple(dataRef.dataId[k] for k in sorted(dataRef.dataId.keys()))
        inputs = dict((key(select.dataRef), select) for select in selectDataList)
        skyMap = patchRef.get(self.config.coaddName + "Coadd_skyMap")
        tract = skyMap[patchRef.dataId["tract"]]
        patch = tract[(tuple(int(i) for i in patchRef.dataId["patch"].split(",")))]
        bbox = patch.getOuterBBox()
        wcs = tract.getWcs()
        cornerPosList = afwGeom.Box2D(bbox).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        dataRefList = self.select.runDataRef(patchRef, coordList, selectDataList=selectDataList).dataRefList
        return [inputs[key(dataRef)] for dataRef in dataRefList]

    def process(self, patchRef, coadd):
        self.log.info("%s: Start processing %s" % (thisNode(), patchRef.dataId))
        try:
            self.processCoadd.process(patchRef, coadd)
        except LsstCppException as e:
            if (isinstance(e.message, InvalidParameterException) and
                re.search("St. dev. must be > 0:", e.message.what())):
                # All the good pixels are outside the area of interest
                self.log.warn("No usable area for detection: %s" % patchRef.dataId)
                return
            raise
        self.log.info("%s: Finished processing %s" % (thisNode(), patchRef.dataId))

    def writeMetadata(self, dataRef):
        pass

"""
StackLauncher:
* Inputs: coadd name (for SkyMap retrieval), PBS stuff, exposure list, reference exposure
* Calculate overlaps: for each input, determine overlapping tract,patch: {visit/ccd --> tract/patch}
* Invert overlaps: {tract/patch --> visit/ccd}
* Determine background reference exposures (for now, use whatever's provided in refId; this is OK for single fields, but eventually we need a means of selecting reference exposures over an extended area; some ideas at https://dev.lsstcorp.org/trac/ticket/2741)
* Persist overlaps and background reference determinations
 save {tract/patch --> ref visit/ccd}
* For each tract/patch, launch StackTask

StackTask:
* Inputs: coadd name, tract/patch, specific exposure list, single ref exposure
* For each input: makeCoaddTempExp
* assembleCoadd with background matching to ref exposure
* (Could do a detection/photometry on the coadd right now, but background subtraction not the best)
* Generate summary statistics required for BackgroundSubtractionTask

BackgroundSubtractionTask:
* Inputs: coadd name, tract, patch list, ref exposure
* Determine background over the patches of the tract that share a common ref exposure id
* For each patch: write background model

StackFit = processCoadd:
* Subtract background, detect and measure
* Measure same apertures for different filters

ResolveDupes = done automatically in forcedPhotCoadd, but processCoadd dupe resolution requires LSST feature (#2893):
* Inputs: All sources measured from coadds
* Resolve duplicates (patch overlaps, tract overlaps)

ForcedPhotometry = forcedPhotCoadd:
* Inputs: coadd name, --id
* Get reference source list, measure with forced apertures

MultiFit = coming....:
* Inputs: coadd name, tract/patch, exposure list
* Get reference source list
* For each source in patch: determine relevant exposures, measure simultaneously across exposures


"""
