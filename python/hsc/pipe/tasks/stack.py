from lsst.pex.config import Config, ConfigurableField
from lsst.pipe.tasks.selectImages import BaseSelectImagesTask, WcsSelectImagesTask
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.base import CmdLineTask


class MpiWcsSelectImagesTask(WcsSelectImagesTask):
    """The selectDataList is only defined properly on the root node, so we need to retrieve it"""
    def runDataRef(self, patchRef, coordList, makeDataRefList=True, selectDataList=[]):
        import pbasf2
        selectDataList = pbasf2.Broadcast(getComm(), selectDataList)
        return super(MpiWcsSelectImagesTask, self).runDataRef(patchRef, coordList, makeDataRefList,
                                                              selectDataList)

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
        return pipeBase.Struct(
            dataRefList = [s.dataRef for s in selectDataList],
            exposureInfoList = [BaseExposureInfo(s.dataRef.dataId, None) for s in selectDataList],
            )


class StackConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name for coadd")
    selectImages = ConfigurableField(target=BaseSelectImagesTask, doc="Select images to process")
    makeCoaddTempExp = ConfigurableField(target=MakeCoaddTempExpTask, doc="Warp images to sky")
    assembleCoadd = ConfigurableField(target=AssembleCoaddTask, doc="Assemble warps into coadd")
    doOverwriteCoadd = Field(dtype=bool, default=False, doc="Overwrite coadd?")

    def setDefaults(self):
        self.selectImages.retarget(MpiWcsSelectImagesTask)
        self.makeCoaddTempExp.selectImages.retarget(NullSelectImagesTask)
        self.assembleCoadd.selectImages.retarget(NullSelectImagesTask)
        self.makeCoaddTempExp.doOverwrite = False

    def validate(self):
        if self.makeCoaddTempExp.coaddName != self.coaddName:
            raise RuntimeError("makeCoaddTempExp.coaddName and coaddName don't match")
        if self.assembleCoadd.coaddName != self.coaddName:
            raise RuntimeError("assembleCoadd.coaddName and coaddName don't match")

class StackTask(PbsCmdLineTask, MpiTask):
    ConfigClass = StackConfig
    _DefaultName = "stack"
    RunnerClass = MpiMultiplexTaskRunner

    def __init__(self, *args, **kwargs):
        super(StackTask, self).__init__(*args, **kwargs)
        self.makeSubtask("selectImages")
        self.makeSubtask("makeCoaddTempExp")
        self.makeSubtask("assembleCoadd")

    @classmethod
    def _makeArgumentParser(cls, doPbs=False):
        """
        Patch references are cheap, so are defined on all nodes.
        Selection references are not cheap (reads Wcs), so are generated
        only on the root node (and only if we're not doing a PBS submission),
        and distributed via the MpiWcsSelectImagesTask.
        """
        parser = MpiArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer, rootOnly=False)
        # We don't want to be reading all the WCSes if we're only in the act of submitting to PBS
        SelectContainerClass = DataIdContainer if doPbs else SelectDataIdContainer
        parser.add_id_argument("--selectId", "raw", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
                               ContainerClass=SelectContainerClass, rootOnly=True)
        return parser

    def run(self, patchRef, selectDataList=[], refDataRef=None):
        """Run the stacking for a single patch

        @param patchRef: Data reference for stack (i.e., contains tract,patch,filter)
        @param selectDataList: List of SelectStruct for inputs
        @param refDataRef: Data reference for reference exposure; None for no background matching
        """
        selectDataList = self.selectExposures(selectDataList)
        self.makeCoaddTempExp.run(patchRef, selectDataList)
        if self.config.doOverwriteCoadd or not patchRef.datasetExists(self.getCoaddDatasetName()):
            coadd = self.assembleCoadd.run(patchRef, selectDataList)
        # XXX Generate and persist summary statistics required for subtracting background over multiple patches
        # XXX Do a basic detection/measurement here? Background subtraction's not the best, but might be useful

    def selectExposures(self, patchRef, selectDataList):
        """Select exposures to operate upon, via the SelectImagesTask

        This is very similar to CoaddBaseTask.selectExposures, except we return
        a list of SelectStruct (same as the input), so we can plug the results into
        future uses of SelectImagesTask.
        """
        key = lambda dataRef: tuple(dataRef.dataId[k] for k in sorted(dataRef.dataId.keys()))
        inputs = dict((key(select.dataRef), select) for select in selectDataList)
        skyInfo = self.getSkyInfo(patchRef)
        cornerPosList = afwGeom.Box2D(skyInfo.bbox).getCorners()
        coordList = [skyInfo.wcs.pixelToSky(pos) for pos in cornerPosList]
        dataRefList = self.select.runDataRef(patchRef, coordList, selectDataList=selectDataList).dataRefList
        return [inputs[key(dataRef)] for dataRef in dataRefList]


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
