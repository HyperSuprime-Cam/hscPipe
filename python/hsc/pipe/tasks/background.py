import argparse
import itertools

import numpy

from lsst.geom import convexHull
from lsst.pex.config import Config, ConfigurableField, ListField, Field
from lsst.pex.exceptions import LsstCppException, DomainErrorException, RuntimeErrorException
from lsst.pipe.base import Task, CmdLineTask, Struct, ArgumentParser
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner, CoaddDataIdContainer, SelectDataIdContainer, getSkyInfo
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
from lsst.pipe.tasks.makeCoaddTempExp import MakeCoaddTempExpTask
from lsst.pipe.tasks.matchBackgrounds import MatchBackgroundsTask
from lsst.pipe.tasks.scaleZeroPoint import ScaleZeroPointTask
from lsst.afw.fits.fitsLib import FitsError
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.afw.cameraGeom.utils as afwCgu

from hsc.pipe.base.mpi import thisNode
from hsc.pipe.base.pool import Pool

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
    return Struct(metadata=md,
                  wcs=wcs,
                  box=box,
                  poly=convexHull([coord.getVector() for coord in imageCorners]),
                  )

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
        return dict((getDataName(dataRef), Struct(dataRef=dataRef,
                                                  poly=imagePoly(dataRef, log=self.log).poly,
                                                  **self.getExposureData(dataRef).getDict())) for
                    dataRef in dataRefList)

    def getExposureData(self, dataRef):
        """Retrieve data about an exposure"""
        self.log.info("Reading exposure data for %s" % (dataRef.dataId,))
        psf = dataRef.get("psf", immediate=True) # requires special handling in Mapper these days
        bgArray = dataRef.get("calexpBackground", immediate=True).getImage().getArray()
        return Struct(psfWidth=psf.computeShape().getTraceRadius(),
                      bgMean=bgArray.mean(),
                      bgStd=bgArray.std(),
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
                visitData = patchOverlaps.get(visitName, None)
                if visitData is None:
                    continue
                selected.append([expData.dataRef for expData in visitData])
            if len(selected) > 0:
                assignments[patch] = selected
        return assignments


class AllInclusiveAssignTask(AssignTask):
    """Select everything we're given.

    A careful selection of a subset may not be necessary when using
    the taper.
    """
    def getExposureData(self, dataRef):
        """No need to read anything"""
        return Struct()
    def calculateBestSelections(self, tractInfo, overlaps, expDataDict):
        """Select everything"""
        return list(set(tuple(expData.dataRef.dataId[k] for k in self.config.visitKeys) for
                        expData in expDataDict.itervalues()))


class TractDataIdContainer(CoaddDataIdContainer):
    """Plug-in to the ArgumentParser to produce a
    list of lists of patch references for each tract.
    """
    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList"""
        coaddName = namespace.config.coaddName + "Coadd"
        datasetType = coaddName + "_bgRef"
        validKeys = namespace.butler.getKeys(datasetType=datasetType, level=self.level)

        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract", "patch"):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(None, "--id must include " + key)

            if "patch" in dataId:
                raise RuntimeError("May not specify 'patch'")

            # tract is required; iterate over it if not provided
            if not "tract" in dataId:
                addList = [dict(tract=tract.getId(), patch="%d,%d" % patch.getIndex(), **dataId)
                           for tract in self.getSkymap(namespace, coaddName) for patch in tract]
            else:
                tract = self.getSkymap(namespace, coaddName)[dataId["tract"]]
                addList = [dict(patch="%d,%d" % patch.getIndex(), **dataId) for patch in tract]

            # Note: adding a list of lists, so we can operate over a single tract
            self.refList += [[namespace.butler.dataRef(datasetType=datasetType, dataId=addId)
                              for addId in addList]]


class CoordinateTransformer(object):
    """A transformation between patch and CCD coordinates

    Allows transformation between the following systems:
    * Patch pixel coordinates (x,y): "patch"
    * Sky coordinates (RA,Dec): "sky"
    * CCD pixel coordinates for one of the CCDs overlapping the patch (i,x,y): "calexp"
    * Physical location on the focal plane (u,v): "position"
    """
    def __init__(self, patchWcs, calexpRefList):
        """Constructor

        @param patchWcs: Wcs for patch
        @param calexpRefList: List of data references for calexps
        """
        self.patchWcs = patchWcs
        self.polyStructList = [imagePoly(ref) for ref in calexpRefList] 
        self.polyList = [each.poly for each in self.polyStructList]
        self.wcsList = [each.wcs for each in self.polyStructList]
        butler = calexpRefList[0].getButler()
        camera = butler.get("camera")
        self.ccdList = [afwCgu.findCcd(camera,
                                       afwCameraGeom.Id(butler.mapper._extractDetectorName(ref.dataId))) for
                        ref in calexpRefList]

    def getImageIndex(self, coord):
        """Return index of ccd that contains provided sky coordinates"""
        vector = coord.getVector()
        for i, poly in enumerate(self.polyList):
            if poly.containsPoint(vector):
                return i
        raise LookupError("No image contains coordinate %s" % coord)

    def skyToCalexp(self, coord):
        index = self.getImageIndex(coord)
        return index, self.wcsList[index].skyToPixel(coord)

    def calexpToSky(self, index, pixel):
        return self.wcsList[index].pixelToSky(pixel)

    def skyToPatch(self, coord):
        return self.patchWcs.skyToPixel(coord)

    def patchToSky(self, pixel):
        return self.patchWcs.pixelToSky(pixel)

    def patchToCalexp(self, pixel):
        return self.skyToCalexp(self.patchToSky(pixel))

    def calexpToPatch(self, index, pixel):
        return self.skyToPatch(self.calexpToSky(index, pixel))

    def patchToPosition(self, pixel):
        try:
            index, ccdPixel = self.patchToCalexp(pixel)
        except LookupError:
            return None
        return self.ccdList[index].getPositionFromPixel(ccdPixel).getMm()

    def getBoresight(self):
        """Calculate the exposure boresight position in sky coordinates

        We average the boresight position provided by each of the detectors.
        This is probably not strictly necessary, but it doesn't cost much.
        """
        coordList = [self.calexpToSky(index, ccd.getPixelFromPosition(afwCameraGeom.FpPoint(0,0))) for
                     index, ccd in enumerate(self.ccdList)]
        boresight = sum((coord.getVector() for coord in coordList), numpy.zeros(3))/len(self.ccdList)
        return afwCoord.Coord(afwGeom.Point3D(boresight))


def cosBell(xx, start, stop):
    """Construct a cos bell, starting and stopping at the nominated positions"""
    ff = numpy.cos(0.5*(xx - start)/(stop - start)*numpy.pi)
    ff[:] = numpy.where(xx < start, 1.0, ff)
    ff[:] = numpy.where(xx > stop, 0.0, ff)
    return ff


class RadiusTaperConfig(Config):
    start = Field(dtype=float, doc="Radius for taper to start")
    stop = Field(dtype=float, doc="Radius for taper to stop")
    def validate(self):
        if self.start < 0 or self.start > self.stop:
            raise RuntimeError("Bad taper radii: %f %f" % (self.start, self.stop))


class RadiusTaperWeightImage(object):
    """Weight image constructed by tapering in radius

    This is a lightweight (low-memory) representation of the
    full weight image, which can be retrieved by calling
    'getImage'.

    Uses a cos bell taper in radius, starting at 'start'
    and cutting off at 'stop'.
    """
    ConfigClass = RadiusTaperConfig

    def __init__(self, bbox, boresight, config):
        """Constructor

        @param bbox: Patch (outer) bounding box
        @param boresight: Boresight position in patch coordinates
        @param config: Configuration
        """
        self._bbox = bbox
        self._boresight = boresight
        self._config = config
        self._taperStart = config.start
        self._taperStop = config.stop

    def getImage(self):
        """Return weight image"""
        image = self.getRadiusImage()
        return self.applyTaper(image)

    def getRadiusImage(self):
        """Return image of radius as a function of position"""
        bbox = self._bbox
        image = afwImage.ImageF(bbox)
        x, y = numpy.ogrid[bbox.getMinY():bbox.getMaxY()+1, bbox.getMinX():bbox.getMaxX()+1]
        xCenter, yCenter = self._boresight
        array = image.getArray()
        array[:] = numpy.sqrt((x - xCenter)**2 + (y - yCenter)**2)
        return image

    def applyTaper(self, image):
        """Apply taper to radius image, in-place"""
        array = image.getArray()
        array[:] = cosBell(array, self._taperStart, self._taperStop)
        return image

class XyTaperConfig(Config):
    xStart = Field(dtype=float, doc="Distance from center for taper in x to start")
    xStop = Field(dtype=float, doc="Distance from center for taper in x to stop")
    yStart = Field(dtype=float, doc="Distance from center for taper in y to start")
    yStop = Field(dtype=float, doc="Distance from center for taper in y to stop")
    def validate(self):
        if self.xStart < 0 or self.xStart > self.xStop:
            raise RuntimeError("Bad taper x distances: %f %f" % (self.xStart, self.xStop))
        if self.yStart < 0 or self.yStart > self.yStop:
            raise RuntimeError("Bad taper y distances: %f %f" % (self.yStart, self.yStop))

class XyTaperWeightImage(object):
    """Weight image constructed by tapering in x and y

    This is a lightweight (low-memory) representation of the
    full weight image, which can be retrieved by calling
    'getImage'.

    Uses a cos bell taper in x and y, each starting and
    cutting off at nominated distances from the boresight.

    Of course, this assumes that the exposure is oriented
    in x,y on the patch.
    """
    ConfigClass = XyTaperConfig
    def __init__(self, bbox, boresight, config):
        """Constructor

        @param bbox: Patch (outer) bounding box
        @param boresight: Boresight position in patch coordinates
        @param config: Configuration
        """
        self._bbox = bbox
        self._boresight = boresight
        self._xStart = config.xStart
        self._xStop = config.xStop
        self._yStart = config.yStart
        self._yStop = config.yStop

    def getImage(self):
        """Return weight image"""
        bbox = self._bbox
        image = afwImage.ImageF(bbox)
        array = image.getArray()
        y, x = numpy.ogrid[bbox.getMinY():bbox.getMaxY()+1, bbox.getMinX():bbox.getMaxX()+1]
        xCenter, yCenter = self._boresight
        array[:] = (cosBell(numpy.abs(x - xCenter), self._xStart, self._xStop)*
                    cosBell(numpy.abs(y - yCenter), self._yStart, self._yStop))
        return image

class ConstructionConfig(Config):
    doWarp = Field(dtype=bool, default=True, doc="Warp images to sky tract/patch?")
    taper = ConfigurableField(target=RadiusTaperWeightImage, doc="Object to provide tapered weight image")
    makeCoaddTempExp = ConfigurableField(target=MakeCoaddTempExpTask, doc="Warp images to sky")
    scaling = ConfigurableField(target=ScaleZeroPointTask, doc="Scale warps to common zero point")
    mask = ListField(dtype=str, default=["BAD", "SAT", "EDGE",], doc="Mask planes to mask")
    matchBackgrounds = ConfigurableField(target=MatchBackgroundsTask, doc="Background matching")
    clobber = Field(dtype=bool, default=False, doc="Clobber existing outputs?")

class ConstructionTask(Task):
    ConfigClass = ConstructionConfig

    def __init__(self, *args, **kwargs):
        super(ConstructionTask, self).__init__(*args, **kwargs)
        self.makeSubtask("makeCoaddTempExp")
        self.makeSubtask("scaling")
        self.makeSubtask("matchBackgrounds")

    def run(self, patchRefList, visitList, overlaps, coaddName="deep"):
        """Construct and write a background reference image from multiple inputs

        This is a persistence layer over the main construction
        method, 'construct'.

        @param patchRefList: List of patch data references
        @param assignmentsList: Embedded lists: patches holding visits holding CCD data references
        @param coaddName: name of coadd
        @return background reference
        """
        pool = Pool()
        datasetType = coaddName + "Coadd_bgRef"
        patchIdList = [patchRef.dataId for patchRef in patchRefList]

        # Check for existing data
        if not self.config.clobber:
            self.log.info("Checking for existing %s data" % datasetType)
            # XXX does this need to be parallel?
            pool.scatterGather(self.checkExisting, True, patchIdList, datasetType)

        # Start with the first visit
        visit = visitList.pop(0)
        pool.scatterGather(self.firstVisit, True, patchIdList, visit=visit, overlaps=overlaps,
                           coaddName=coaddName)

        # Add in each visit one by one
        for visit in visitList:
            # Generate background matching model
            bgModelList = pool.scatterGatherToPrevious(self.matchNextVisit, patchIdList, visit=visit,
                                                       overlaps=overlaps, coaddName=coaddName)

            # Merge background models
            bgModel = self.mergeBackgroundModels(bgModelList)

            # Use common background model when adding the next visit
            pool.scatterGatherToPrevious(self.addNextVisit, patchIdList, visit=visit, overlaps=overlaps,
                                         cooadName=coaddName)

        # Finish up and write
        pool.scatterGatherToPrevious(self.finalize, patchIdList, datasetType=datasetType)

    def checkExisting(self, cache, patchId, datasetType):
        butler = cache.butler
        patchRef = getDataRef(butler, patchId)
        if patchRef.datasetExists(datasetType):
            raise RuntimeError("Background reference exists for %s" % (patchRef.dataId,))

    def firstVisit(self, cache, patchId, visit, overlaps, coaddName):
        patchRef = getDataRef(cache.butler, patchId)
        XXXX

    def construct(self, patchRef, assignments, coaddName="deep"):
        """Construct a background reference image from multiple inputs

        @param patchRef: data reference for a patch
        @param assignments: list of a list of data references for each visit
        @param coaddName: name of coadd
        @return background reference
        """
        warpRefList = []
        weightList = []
        skyInfo = getSkyInfo(coaddName, patchRef)
        for calexpRefList in assignments:
            warpRef = self.getWarpRef(patchRef, calexpRefList, coaddName=coaddName)
            try:
                if self.config.doWarp:
                    self.warp(warpRef, calexpRefList, skyInfo, coaddName=coaddName)
                weight = self.generateWeight(skyInfo, calexpRefList)
            except RuntimeError as e:
                self.log.warn("Ignoring warp %s due to %s: %s" % (warpRef.dataId, e.__class__.__name__, e))
                continue
            warpRefList.append(warpRef)
            weightList.append(weight)

        num = len(warpRefList)
        if num == 0:
            raise RuntimeError("No good input exposures")

        bg = None
        for warpRef, weight in zip(warpRefList, weightList):
            self.log.info("Adding warp %s" % (warpRef.dataId,))
            bg = self.addWarp(bg, warpRef, weight, coaddName=coaddName)

        bgImage = bg.exposure.getMaskedImage()
        bgImage /= bg.weight
        return bg.exposure

    def warp(self, warpRef, calexpRefList, skyInfo, coaddName="deep"):
        """Warp CCDs to patch

        @param warpRef: data reference for warp
        @param calexpRefList: list of calexp data references
        @param skyInfo: struct with skymap information
        @param coaddName: name of coadd
        """
        datasetType = coaddName + "Coadd_tempExp"
        if warpRef.datasetExists(datasetType):
            return
        exp = self.makeCoaddTempExp.createTempExp(calexpRefList, skyInfo)
        warpRef.put(exp, datasetType)

    def generateWeight(self, skyInfo, calexpRefList):
        """Construct a weight

        Returns a lightweight version of the actual weight image;
        call 'getImage' for the actual image.

        @param skyInfo: struct with skymap information
        @param calexpRefList: list of calexp data references
        @return light weight image
        """
        transformer = CoordinateTransformer(skyInfo.wcs, calexpRefList)
        boresight = transformer.getBoresight()
        return self.config.taper.apply(skyInfo.bbox, transformer.skyToPatch(boresight))

    def getWarpRef(self, patchRef, calexpRefList, coaddName="deep"):
        """Generate a warp data reference

        A "warp" is also known as a "coaddTempExp".

        @param patchRef: data reference for the patch
        @param calexpRefList: list of data references for constituent calexps
        @return warp data reference
        """
        dataId = patchRef.dataId.copy()
        dataId.update(calexpRefList[0].dataId)
        for calexpRef in calexpRefList[1:]:
            for key, value in calexpRef.dataId.iteritems():
                if key in dataId and dataId[key] != value:
                    del dataId[key]

        tempExpName = coaddName + "Coadd_tempExp"
        return patchRef.getButler().dataRef(datasetType=tempExpName, dataId=dataId)

    def matchWarp(self, bg, warpRef, weight, coaddName="deep"):
        """Match a new warp to the background reference

        @param bg: background reference struct (exposure,weight elements)
        @param warpRef: warp data reference
        @param weight: weight image (or lightweight version)
        @return background model
        """
        warp = warpRef.get(coaddName + "Coadd_tempExp", immediate=True)
        self.scaling.computeImageScaler(warp, warpRef).scaleMaskedImage(warp.getMaskedImage())
        if hasattr(weight, "getImage"):
            weight = getattr(weight, "getImage")()

        # Remove masked pixels
        mask = warp.getMaskedImage().getMask()
        bitmask = mask.getPlaneBitMask(self.config.mask)
        weightArray = weight.getArray()
        warpArray = warp.getMaskedImage().getImage().getArray()
        varArray = warp.getMaskedImage().getVariance().getArray()
        bad = numpy.logical_or(mask.getArray() & bitmask > 0,
                               numpy.logical_not(numpy.isfinite(warpArray)))
        warpArray[bad] = 0.0
        weightArray[bad] = 0.0
        varArray[bad] = 0.0
        del bad

        # Chop off overlaps so we don't step on other patches' toes
        bbox = patchInfo.getInnerBBox()
        bgImage = afwImage.MaskedImageF(bg.exposure.getMaskedImage(), bbox).clone()
        subWeight = afwImage.ImageF(bg.weight, bbox)

        # Undo weighting of background, so we can match
        bgImage /= subWeight

        bgImage.__isub__(warp.getMaskedImage())
        bgModel = afwMath.BackgroundMIF(bgImage, bgCtrl)
        return bgModel



class BackgroundReferenceConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd of interest")
    select = ConfigurableField(target=WcsSelectImagesTask, doc="Task to select input images")
    assign = ConfigurableField(target=AssignTask, doc="Task to assign inputs to patches")
    construct = ConfigurableField(target=ConstructionTask, doc="Task to construct background reference")

class BackgroundReferenceTask(CmdLineTask):
    _DefaultName = "bgRef"
    ConfigClass = BackgroundReferenceConfig
    RunnerClass = CoaddTaskRunner

    def __init__(self, config, **kwargs):
        super(BackgroundReferenceTask, self).__init__(config=config, **kwargs)
        self.makeSubtask("select")
        self.makeSubtask("assign")
        self.makeSubtask("construct")
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

    def run(self, patchRefList, selectDataList=[]):
        skyMap = patchRefList[0].get(self.config.coaddName + "Coadd_skyMap")
        tract = skyMap[patchRefList[0].dataId['tract']]
        dataRefList = self.selectExposures(tract, selectDataList)
        assignments = self.assign.run(tract, dataRefList)
        self.constructBackgrounds(patchRefList, assignments)

    def constructBackgrounds(self, patchRefList, assignments):
        """Construct background reference images

        This method could be trivially parallelised.

        @param patchRefList: List of patch data references
        @param assignments: List of a list of data references for each visit
        """
        for patchRef in patchRefList:
            patch = tuple(map(int, patchRef.dataId["patch"].split(",")))
            try:
                self.construct.run(patchRef, assignments[patch], coaddName=self.config.coaddName)
            except RuntimeError as e:
                self.log.warn("Unable to construct background reference for %s due to %s: %s" %
                              (patch, e.__class__.__name__, e))
                continue

    def selectExposures(self, tractInfo, selectDataList=[]):
        """Select exposures to include

        @param tractInfo: tract object from skymap
        @param selectDataList: list of SelectStruct
        @return list of calexp data references
        """
        wcs = tractInfo.getWcs()
        cornerPosList = afwGeom.Box2D(tractInfo.getBBox()).getCorners()
        coordList = [wcs.pixelToSky(pos) for pos in cornerPosList]
        return self.select.runDataRef(None, coordList, selectDataList=selectDataList).dataRefList

    def writeMetadata(self, *args, **kwargs):
        pass
    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchemas(self, *args, **kwargs):
        pass


class MpiBackgroundReferenceConfig(BackgroundReferenceConfig):
    """We're using this as part of the super-stacker, so the warps have already
    been constructed and we don't want to risk regenerating them with a different
    configuration.
    """
    makeCoaddTempExp = None
    def setDefaults(self):
        super(MpiBackgroundReferenceConfig, self).setDefaults()
        self.construct.doWarp = False

class MpiBackgroundReferenceTask(MpiTask, BackgroundReferenceTask):
    """MPI-enabled background reference construction, intended for use within the super-stacker"""
    ConfigClass = MpiBackgroundReferenceConfig

    def run(self, patchRefList, selectDataList=[]):
        """Construct a set of background references

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param patchRefList: List of patch data references
        @param selectDataList: List of SelectStruct for inputs
        """
        if self.rank == self.root:
            self.log.info("%s: Root node calculating overlaps and assignment" % thisNode())
            super(MpiBackgroundReferenceTask, self).run(patchRefList, selectDataList=selectDataList)
        else:
            # Must get the slave nodes into the same place the master node ends up
            self.constructBackgrounds(patchRefList, None)

    def constructBackgrounds(self, patchRefList, assignments):
        """Farm out the construction of background references

        All nodes execute this method, though the master and slaves
        take different routes through it.

        @param patchRefList: List of patch data references
        @param assignments: List of a list of data references for each visit
        """
        import pbasf2
        if self.rank == self.root:
            args = [Struct(patchRef=patchRef,
                           assignment=assignments.get(tuple(map(int, patchRef.dataId["patch"].split(","))),
                                                      None),
                           ) for patchRef in patchRefList]
        else:
            args = None
        self.log.info("%s: Ready to construct backgrounds" % thisNode())
        pbasf2.ScatterJob(self.comm, self.constructSingleBackground, args, root=self.root)

    def constructSingleBackground(self, struct):
        """Wrapper for ConstructBackgroundTask

        Only slave nodes execute this method.

        Because only one argument may be passed, it is expected to
        contain multiple elements, which are:
        @param patchRef: data reference for patch
        @param assignments: List of data references for each visit
        """
        patchRef = struct.patchRef
        assignment = struct.assignment
        self.log.info("%s: Start constructing background for %s" % (thisNode(), patchRef.dataId))
        if assignment is None:
            self.log.warn("%s: No inputs assigned for %s" % (thisNode(), patchRef.dataId))
            return
        try:
            self.construct.run(patchRef, assignment, coaddName=self.config.coaddName)
        except RuntimeError as e:
            self.log.warn("Unable to construct background reference for %s due to %s: %s" %
                          (patchRef.dataId, e.__class__.__name__, e))
        self.log.info("%s: Finished constructing background for %s" % (thisNode(), patchRef.dataId))


"""
This doesn't work as I've approached it: the matching is being performed patch by patch instead of
consistently over the entire tract.  For example, an edge patch that doesn't have the same first warp as the
others will end up matching to a different warp, resulting in a discontinuity across the patch boundary.  More
subtle discontinuities will result from slightly different and discontinuous background matching solutions for
neighbouring patches.

What is required is to do the background matching solution over the entire tract at once.  The subtraction
should be performed patch by patch on the slaves, and the background-difference samples returned to the master
node for a single, continuous background-difference model to be constructed and returned to the slaves for
application.  This is likely going to require some extra hooks in the background model code to merge multiple
models from different parts of a larger image, and to allow pickling.
"""
