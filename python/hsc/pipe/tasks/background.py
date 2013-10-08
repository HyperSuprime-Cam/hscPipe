import os
import argparse
import itertools
from collections import Counter
from lsst.geom import convexHull
from lsst.pex.config import Config, ConfigurableField, ListField, Field
from lsst.pex.exceptions import LsstCppException, DomainErrorException, RuntimeErrorException
from lsst.pipe.base import Task, CmdLineTask, Struct, ArgumentParser
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner, CoaddDataIdContainer, SelectDataIdContainer
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
from lsst.afw.fits.fitsLib import FitsError
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage


def imagePoly(dataRef, log=None, imageName="calexp"):
    """Generate a SphericalConvexPolygon from an image

    If a polygon cannot be generated, returns None.

    @param dataRef: data reference
    @param log: logger, for errors
    @param imageName: butler data type to retrieve
    @return SphericalConvexPolygon for image
    """
    try:
        md = dataRef.get(imageName + "_md", immediate=True)
    except LsstCppException as e:
        if not isinstance(e.message, FitsError): # Unable to open file
            raise
        if log: log.warn("Unable to read calexp header for %s" % (dataRef.dataId))
        return None

    wcs = afwImage.makeWcs(md)
    box = afwGeom.Box2D(afwGeom.Point2D(0,0), afwGeom.Extent2D(md.get("NAXIS1"), md.get("NAXIS2")))
    try:
        imageCorners = [wcs.pixelToSky(pix) for pix in box.getCorners()]
    except LsstCppException as e:
        # Protecting ourselves from awful Wcs solutions in input images
        if (not isinstance(e.message, DomainErrorException) and
            not isinstance(e.message, RuntimeErrorException)):
            raise
        if log: log.logdebug("WCS error in testing calexp %s: %s" % (dataRef.dataId, e))
        return None
    return convexHull([coord.getVector() for coord in imageCorners])

def getDataName(dataRef):
    """Generate a hashable name for a dataRef

    This implementation assumes that all dataRefs have the
    same set of keys.
    """
    return tuple(dataRef.dataId[k] for k in sorted(dataRef.dataId.keys()))

class OverlappingVisit(Struct):
    """Container for a visit overlapping a patch

    A visit is a group of exposures (typically the CCDs comprising the visit)
    that overlap a common patch.

    The overlap fraction is accumulated for all component exposures.
    """
    def __init__(self, visit):
        super(OverlappingVisit, self).__init__(visit=visit, _components=[], overlap=0.0)

    def add(self, expOverlap):
        """Add an overlapping exposure"""
        assert expOverlap.visit == self.visit
        self._components.append(expOverlap)
        self.overlap += expOverlap.overlap

    def __iter__(self):
        return iter(self._components)

class OverlappingExposure(Struct):
    """Container for an exposure (typically a CCD) overlapping a patch"""
    def __init__(self, dataRef, visit, intersection, patchArea):
        """Constructor

        @param dataRef: data reference
        @param visit: Name (identifying tuple) of visit
        @param itersection: SphericalConvexPolygon of intersection between exposure and patch
        @param patchArea: Area of patch, for calculating the fractional overlap
        """
        super(OverlappingExposure, self).__init__(dataRef=dataRef, dataName=getDataName(dataRef),
                                                  visit=visit, intersection=intersection,
                                                  overlap=intersection.area()/patchArea)

class AssignConfig(Config):
    visitKeys = ListField(dtype=str, default=["visit"], doc="dataId keys that identify an exposure")
    scoreNumExposures = Field(dtype=float, default=2.0, doc="Weight for scoring number of exposures")
    scoreNumPatches = Field(dtype=float, default=2.0, doc="Weight for scoring number of patches")
    scoreAverageOverlap = Field(dtype=float, default=1.0, doc="Weight for scoring average overlap")
    scoreAverageBgStdev = Field(dtype=float, default=0.01, doc="Weight for scoring average background stdev")
    scoreAveragePsfWidth = Field(dtype=float, default=0.1, doc="Weight for scoring average PSF width")

class AssignTask(Task):
    ConfigClass = AssignConfig

    def run(self, tractInfo, dataRefList):
        """Assign a set of exposures to each patch for background reference creation

        @param tractInfo: tract for which to select
        @param dataRefList: list of input (calexp) data references
        @return dict mapping patch (x,y) tuple to a list of a list of data references for each visit
        """
        expDataDict = self.gatherExposureData(dataRefList)
        overlaps = self.calculateOverlaps(tractInfo, expDataDict)
        selections = self.calculateBestSelections(tractInfo, overlaps, expDataDict)
        assignments = self.extractSelections(selections, overlaps)
        return assignments

    def calculateOverlaps(self, tractInfo, expDataDict):
        """Calculate overlaps between all patches and all input exposures

        The returned overlaps is a dict mapping patches (by x,y tuple) to
        a dict of overlapping visits.  The dict of overlapping visits maps
        the visit name (tuple of dataId values) to an OverlappingVisit
        object, which contains a list of OverlappingExposure objects.

        @param tractInfo: tract containing patches
        @param expDataDict: input data (dict mapping data name --> struct with data)
        @return overlaps (see above)
        """
        overlaps = {}
        tractWcs = tractInfo.getWcs()
        xNum, yNum = tractInfo.getNumPatches()
        for y in range(yNum):
            for x in range(xNum):
                patch = tractInfo.getPatchInfo((x, y))
                patchCorners = afwGeom.Box2D(patch.getOuterBBox()).getCorners()
                patchVertices = [tractWcs.pixelToSky(cnr).getVector() for cnr in patchCorners]
                patchPoly = convexHull(patchVertices)
                patchOverlaps = {}
                for expData in expDataDict.itervalues():
                    intersection = patchPoly.intersect(expData.poly)
                    if intersection is None:
                        continue
                    area = intersection.area()
                    visitName = tuple(expData.dataRef.dataId[k] for k in self.config.visitKeys)
                    if visitName in patchOverlaps:
                        visitOverlaps = patchOverlaps[visitName]
                    else:
                        visitOverlaps = OverlappingVisit(visitName)
                        patchOverlaps[visitName] = visitOverlaps
                    visitOverlaps.add(OverlappingExposure(expData.dataRef, visitName, intersection,
                                                          patchPoly.area()))
                if patchOverlaps:
                    overlaps[(x,y)] = patchOverlaps

        return overlaps

    def gatherExposureData(self, dataRefList):
        """Gather data about each input exposure

        This is currently implemented as a serial operation, but could be
        made parallel using MPI.

        @param dataRefList: list of (calexp) data references
        @return dict mapping data name to struct with exposure data
        """
        import pickle
        if False:
            expDataDict = dict((getDataName(dataRef), self.getExposureData(dataRef)) for
                                dataRef in dataRefList)
            pickle.dump(expDataDict, open("expData.pickle", "w"))
        else:
            expDataDict = pickle.load(open("expData.pickle", "r"))

        return expDataDict

    def getExposureData(self, dataRef):
        """Retrieve data about an exposure"""
        self.log.info("Reading exposure data for %s" % (dataRef.dataId,))
        psf = dataRef.get("psf", immediate=True) # requires special handling in Mapper these days
        bgArray = dataRef.get("calexpBackground", immediate=True).getImage().getArray()
        return Struct(dataRef=dataRef,
                      poly=imagePoly(dataRef, log=self.log),
                      psfWidth=psf.computeShape().getTraceRadius(),
                      bgMean=bgArray.mean(),
                      bgStd=bgArray.std()
                      )

    def scoreAssignments(self, assignments, selectedData, selections, overlaps):
        """Provide a score (metric) for a possible set of assignments

        The sense of the score is that a higher score is better.

        This uses a naive metric that will probably need tuning and/or replacement.

        @param assignments: XXX
        @param selectedData: dict mapping of data name to exposure data, only for selected exposures
        @param selections: list of visit names (dataId subset tuples) that have been selected
        @param overlaps: overlaps between patches and exposures
        @return score
        """
        score = 0

        # Want few distinct exposures involved:
        # ==> smaller number of exposures --> higher score
        score += self.config.scoreNumExposures/len(selections)

        # Want exposures to cover largest possible area:
        # ==> larger number of patches per exposure --> higher score
        score += self.config.scoreNumPatches*len(assignments)/len(overlaps)

        # Want patch overlaps to be substantial:
        # ==> larger patch overlaps --> higher score
        averageOverlap = sum(max(visitData.overlap for visitData in patchAssignments.itervalues()) for
                             patchAssignments in assignments.itervalues())/len(assignments)
        score += self.config.scoreAverageOverlap*averageOverlap

        # Want background to be fairly constant:
        # ==> smaller background standard deviation --> higher score
        averageBgStdev = sum(data.bgStd for data in selectedData.itervalues())/len(assignments)
        score += self.config.scoreAverageBgStdev/averageBgStdev

        # Want PSF to be decent (easier subtractions, smaller kernel when PSF-matching):
        # ==> smaller PSF --> higher score
        averagePsfWidth = sum(data.psfWidth for data in selectedData.itervalues())/len(assignments)
        score += self.config.scoreAveragePsfWidth/averagePsfWidth

        return score

    def scoreSelections(self, selections, overlaps, expDataDict):
        """Provide a score (metric) for a possible set of selections

        This method calculates the assignments of exposures to patches
        and then passes on the work of calculating the score to
        scoreAssignments.

        @param selections: list of visit names (dataId subset tuples) that have been selected
        @param overlaps: overlaps between patches and exposures
        @param expDataDict: dict mapping data name to exposure data
        @return score
        """
        assignments = {}
        selectedData = {}
        for patch, patchOverlaps in overlaps.iteritems():
            selected = {}
            for visitName, visitOverlaps in patchOverlaps.iteritems():
                if visitName not in selections:
                    continue
                selected[visitName] = visitOverlaps

                for expOverlaps in visitOverlaps:
                    dataName = expOverlaps.dataName
                    if not dataName in selectedData:
                        selectedData[dataName] = expDataDict[dataName]
            if len(selected) > 0:
                assignments[patch] = selected

        return self.scoreAssignments(assignments, selectedData, selections, overlaps)

    def calculateBestSelections(self, tractInfo, overlaps, expDataDict):
        """Calculate the best selections of visits

        This implementation scores all possible selections and returns
        the best.  This may be suitable for small tracts, but for larger
        tracts (with more combinations) something else may be required.

        @param tractInfo: tract with patches
        @param overlaps: overlaps between patches and exposures
        @param expDataDict: dict mapping data name to exposure data
        @return list of visit names (dataId subset tuples)
        """
        bestScore = 0
        bestSelections = None

        def selectionsIterator(overlaps):
            """A generator function for all possible selections"""
            exposureList = set()
            for patchOverlaps in overlaps.itervalues():
                exposureList.update(set(patchOverlaps.keys()))

            for num in range(len(exposureList)):
                for expList in itertools.combinations(exposureList, num + 1):
                    yield expList

        for selections in selectionsIterator(overlaps):
            score = self.scoreSelections(selections, overlaps, expDataDict)

            self.log.info("Selection %s ==> score = %f" % (selections, score))

            if False:
                self.visualizeSelections(tractInfo, selections, overlaps, expDataDict)
            if score > bestScore:
                bestScore = score
                bestSelections = selections

        self.log.info("Best selection is %s ==> score = %f" % (bestSelections, bestScore))
        return bestSelections

    def visualizeSelections(self, tractInfo, selections, overlaps, expDataDict):
        """Visualise selections using matplotlib

        Plots the patches and the selected exposures overlaid.

        @param tractInfo: tract with patches
        @param selections: list of visit names (dataId subset tuples)
        @param overlaps: overlaps between patches and exposures
        @param expDataDict: dict mapping data name to exposure data
        """
        import matplotlib
        import matplotlib.pyplot as plt

        fig = plt.figure()
        plt.clf()
        axes = fig.add_subplot(1,1,1)
        colorCycle = matplotlib.rcParams['axes.color_cycle']

        wcs = tractInfo.getWcs()

        for patchIndex, patch in enumerate(tractInfo):
            box = patch.getInnerBBox()

            xMin, xMax, yMin, yMax = box.getMinX(), box.getMaxX(), box.getMinY(), box.getMaxY()
            axes.plot((xMin, xMax, xMax, xMin, xMin), (yMin, yMin, yMax, yMax, yMin), 'k:')

        for i, visitName in enumerate(selections):
            color = colorCycle[i % len(colorCycle)]
            for patchOverlaps in overlaps.itervalues():
                if not visitName in patchOverlaps:
                    continue
                for expOverlap in patchOverlaps[visitName]:
                    expData = expDataDict[expOverlap.dataName]
                    vertices = expData.poly.getVertices()
                    corners = [wcs.skyToPixel(afwCoord.Coord(afwGeom.Point3D(x,y,z))) for x,y,z in vertices]
                    xCorners = [cnr.getX() for cnr in corners]
                    yCorners = [cnr.getY() for cnr in corners]
                    axes.plot(xCorners + [xCorners[0]], yCorners + [yCorners[0]], color + '-')

        plt.show()

    def extractSelections(self, selections, overlaps):
        """Extract the selected data references

        @param selections: list of visit names (dataId subset tuples)
        @param overlaps: overlaps between patches and exposures
        @return dict mapping patch x,y tuple to a list lists of exposure data references for each visit
        """
        assignments = {}
        for patch, patchOverlaps in overlaps.iteritems():
            selected = []
            for visitName in selections:
                visitData = patchOverlaps[visitName]
                selected.append([expData.dataRef for expData in visitData])
            if len(selected) > 0:
                assignments[patch] = selected
        return assignments


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
        name = self._filename(tractRef, coaddName)
        dirname = os.path.dirname(name)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        f = open(name, "w")
        self._pickle.dump(assignments, f)

    def read(self, tractRef, coaddName):
        f = open(self._filename(tractRef, coaddName))
        return self._pickle.load(f)


class BackgroundReferenceConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd of interest")
    io = ConfigurableField(target=BackgroundReferenceIoTask, doc="I/O for background references")
    select = ConfigurableField(target=WcsSelectImagesTask, doc="Task to select input images")
    assign = ConfigurableField(target=AssignTask, doc="Task to assign inputs to patches")
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
        parser.add_id_argument("--selectId", "calexp", help="data ID, e.g. --selectId visit=6789 ccd=0..9",
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
            skyMap = tractRef.get(self.config.coaddName + "Coadd_skyMap")
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
