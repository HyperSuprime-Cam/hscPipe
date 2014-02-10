import math
import argparse
import itertools

import numpy

from lsst.geom import convexHull
from lsst.pex.config import Config, ConfigurableField, ListField, Field, ConfigField, ChoiceField
from lsst.pex.exceptions import LsstCppException, DomainErrorException, RuntimeErrorException
from lsst.pipe.base import Task, Struct, ArgumentParser
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner, CoaddDataIdContainer, SelectDataIdContainer, getSkyInfo
from lsst.pipe.tasks.selectImages import WcsSelectImagesTask
from lsst.pipe.tasks.scaleZeroPoint import ScaleZeroPointTask
from lsst.afw.fits.fitsLib import FitsError
import lsst.afw.display.ds9 as ds9
import lsst.afw.math as afwMath
import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCameraGeom
import lsst.afw.cameraGeom.utils as afwCgu
import lsst.meas.algorithms as measAlg

from hsc.pipe.base.pool import Pool, Debugger, NODE
from hsc.pipe.base.pbs import PbsPoolTask

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
    scoreNumExposures = Field(dtype=float, default=1.0, doc="Weight for scoring number of exposures")
    scoreNumPatches = Field(dtype=float, default=2.0, doc="Weight for scoring number of patches")
    scoreAverageOverlap = Field(dtype=float, default=5.0, doc="Weight for scoring average overlap")
    scoreAverageBgStdev = Field(dtype=float, default=0.01, doc="Weight for scoring average background stdev")
    scoreAveragePsfWidth = Field(dtype=float, default=0.1, doc="Weight for scoring average PSF width")


def getPatchIndex(dataId):
    return tuple(map(int, dataId["patch"].split(",")))

class AssignTask(Task):
    ConfigClass = AssignConfig

    def run(self, tractInfo, dataRefList, visitKeys=["visit"]):
        """Assign a set of exposures to each patch for background reference creation

        @param tractInfo: tract for which to select
        @param dataRefList: list of input (calexp) data references
        @return dict mapping patch (x,y) tuple to a list of a list of data references for each visit
        """
        pool = Pool(None)
        expData = pool.map(self.gatherExposureData, dataRefList)
        expDataDict = {}
        for dataRef, data in zip(dataRefList, expData):
            data.dataRef = dataRef
            expDataDict[getDataName(dataRef)] = data
        overlaps = self.calculateOverlaps(tractInfo, expDataDict, visitKeys)
        selections = self.calculateBestSelections(tractInfo, overlaps, expDataDict, visitKeys)
        assignments = self.extractSelections(selections, overlaps)
        return Struct(overlaps=overlaps, visits=selections, assignments=assignments)

    def calculateOverlaps(self, tractInfo, expDataDict, visitKeys=["visit"]):
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
                    visitName = tuple(expData.dataRef.dataId[k] for k in visitKeys)
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

    def gatherExposureData(self, dataRef):
        """Gather data about an input exposure

        This method adds a polygon, which is essential for determining
        which patches the exposure overlaps.  The "getExposureData"
        method adds other data, useful for scoring.

        @param dataRef: data reference (for calexp)
        @return struct with exposure data
        """
        return Struct(poly=imagePoly(dataRef, log=self.log).poly, **self.getExposureData(dataRef).getDict())

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
        # Want few distinct exposures involved:
        # ==> smaller number of exposures --> higher score
        scoreNum = self.config.scoreNumExposures/len(selections)

        # Want exposures to cover largest possible area:
        # ==> larger number of patches per exposure --> higher score
        scoreArea = self.config.scoreNumPatches*len(assignments)/len(overlaps)

        # Want patch overlaps to be substantial:
        # ==> larger patch overlaps --> higher score
        averageOverlap = sum(max(visitData.overlap for visitData in patchAssignments.itervalues()) for
                             patchAssignments in assignments.itervalues())/len(overlaps)
        scoreOverlap = self.config.scoreAverageOverlap*averageOverlap

        # Want background to be fairly constant:
        # ==> smaller background standard deviation --> higher score
        averageBgStdev = sum(data.bgStd for data in selectedData.itervalues())/len(assignments)
        scoreBgStdev = self.config.scoreAverageBgStdev/averageBgStdev

        # Want PSF to be decent (easier subtractions, smaller kernel when PSF-matching):
        # ==> smaller PSF --> higher score
        averagePsfWidth = sum(data.psfWidth for data in selectedData.itervalues())/len(assignments)
        scorePsfWidth = self.config.scoreAveragePsfWidth/averagePsfWidth

        self.log.info("Score %s: %f %f %f %f %f" % (selections, scoreNum, scoreArea, scoreOverlap,
                                                    scoreBgStdev, scorePsfWidth))
        return scoreNum + scoreArea + scoreOverlap + scoreBgStdev + scorePsfWidth

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

    def calculateBestSelections(self, tractInfo, overlaps, expDataDict, visitKeys=["visit"]):
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
        if bestSelections is None:
            raise RuntimeError("Unable to select background reference visits")
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
    def calculateBestSelections(self, tractInfo, overlaps, expDataDict, visitKeys=["visit"]):
        """Select everything"""
        return list(set(tuple(expData.dataRef.dataId[k] for k in visitKeys) for
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


class MatchBackgroundsConfig(Config):
    background = ConfigField(dtype=measAlg.BackgroundConfig, doc="Background matching config")
    undersampleStyle = ChoiceField(
        doc="behaviour if there are too few points in grid for requested interpolation style",
        dtype=str, default="REDUCE_INTERP_ORDER",
        allowed={
            "THROW_EXCEPTION": "throw an exception if there are too few points",
            "REDUCE_INTERP_ORDER": "use an interpolation style with a lower order.",
            "INCREASE_NXNYSAMPLE": "Increase the number of samples used to make the interpolation grid.",
            }
        )
    interpolation = ChoiceField(
        doc="how to interpolate the background values. This maps to an enum; see afw::math::Background",
        dtype=str, default="AKIMA_SPLINE",
        allowed={
            "CONSTANT" : "Use a single constant value",
            "LINEAR" : "Use linear interpolation",
            "NATURAL_SPLINE" : "cubic spline with zero second derivative at endpoints",
            "AKIMA_SPLINE": "higher-level nonlinear spline that is more robust to outliers",
            "NONE": "No background estimation is to be attempted",
            }
        )

class MatchBackgroundsTask(Task):
    ConfigClass = MatchBackgroundsConfig

    def run(self, ref, warp, inPlace=True):
        """Generate a background matching model

        @param ref: Background reference image
        @param warp: Science warp image
        @param inPlace: Operate on background reference in-place?
        @return background model
        """
        diff = self.subtract(ref, warp, inPlace=inPlace)
        bgModel = self.calculateBackgroundModel(diff)
        return bgModel

    def _getImage(self, image):
        """Get an image suitable for arithmetic operations

        Mainly concerned with getting a MaskedImage out of an Exposure.
        """
        if hasattr(image, "getMaskedImage"):
            return image.getMaskedImage()
        return image

    def subtract(self, ref, warp, inPlace=True):
        """Subtract warp from reference

        This is drop-dead simple now; might want to do PSF-matching.

        @param ref: Background reference image
        @param warp: Science warp image
        @param inPlace: Operate on background reference in-place?
        @return Subtracted image
        """
        if not inPlace:
            ref = ref.clone()
        refImage = self._getImage(ref)
        warpImage = self._getImage(warp)
        refImage -= warpImage
        return ref

    def calculateBackgroundModel(self, diff):
        """Determine a background model from a diff image

        This is the heart of background matching.

        Bad pixels are zeroed out to ensure the background matching
        addition goes to zero where there are no pixels.

        @param diff: Difference exposure
        @return background model
        """
        diff = self._getImage(diff)
        bgModel = measAlg.getBackground(diff, self.config.background)
        return afwMath.cast_BackgroundMI(bgModel)

    def apply(self, warp, bgModel, inPlace=True):
        """Apply background matching model to warp

        @param warp: Science warp image
        @param bgModel: Background model to apply
        @param inPlace: Operate on warp in-place?
        @return Background-matched warp
        """
        if not inPlace:
            warp = warp.clone()
        warpImage = self._getImage(warp)
        bgImage = bgModel.getImageF(self.config.interpolation, self.config.undersampleStyle)
        warpImage += bgImage
        return warpImage


class ConstructionConfig(Config):
    taper = ConfigurableField(target=RadiusTaperWeightImage, doc="Object to provide tapered weight image")
    scaling = ConfigurableField(target=ScaleZeroPointTask, doc="Scale warps to common zero point")
    mask = ListField(dtype=str, default=["BAD", "SAT", "EDGE",], doc="Mask planes to mask")
    matching = ConfigurableField(target=MatchBackgroundsTask, doc="Background matching")
    bgSize = Field(dtype=int, default=256, doc="Merged background model bin size (tract pixels)")
    clobber = Field(dtype=bool, default=False, doc="Clobber existing outputs?")

class ConstructionTask(Task):
    ConfigClass = ConstructionConfig

    def __init__(self, *args, **kwargs):
        super(ConstructionTask, self).__init__(*args, **kwargs)
        self.makeSubtask("scaling")
        self.makeSubtask("matching")
        self.debug = False

    def run(self, patchRefList, visitList, overlaps, coaddName="deep", visitKeys=["visit"]):
        """Construct and write a background reference image from multiple inputs

        This is a persistence layer over the main construction
        method, 'construct'.

        @param patchRefList: List of patch data references
        @param visitList: List of visit names
        @param overlaps: Dict of patch overlaps
        @param coaddName: name of coadd
        @param visitKeys: Name of keys that define a visit
        """
        warpType = coaddName + "Coadd_tempExp"
        bgRefType = coaddName + "Coadd_bgRef"

        butler = patchRefList[0].getButler()

        pool = Pool("background")
        pool.storeSet(butler=butler, coaddName=coaddName, visitKeys=visitKeys,
                      warpType=warpType, bgRefType=bgRefType)

        # Check for existing data
        if not self.config.clobber:
            self.log.info("Checking for existing %s data" % bgRefType)
            for patchRef in patchRefList:
                if patchRef.datasetExists(bgRefType):
                    raise RuntimeError("Background reference exists for %s" % (patchRef.dataId,))


        def extractPatchData(visit):
            patchDataList = []
            for patchRef in patchRefList:
                patchIndex = getPatchIndex(patchRef.dataId)
                patchOverlaps = overlaps[patchIndex]
                calexpRefList = ([comp.dataRef for comp in patchOverlaps[visit]] if
                                 visit in patchOverlaps else None)
                patchDataList.append(Struct(patchRef=patchRef, calexpRefList=calexpRefList))
            return patchDataList

        # Start with the first visit
        patchIdList = [patchRef.dataId for patchRef in patchRefList]
        visitList = list(visitList) # Our own, mutable copy
        visit = visitList.pop(0)
        pool.map(self.firstVisit, extractPatchData(visit), visit)

        # Add in each visit one by one, using a common background model
        for visit in visitList:
            bgModelList = pool.mapToPrevious(self.matchNextVisit, patchRefList, visit=visit)
            bgModel = self.mergeBackgroundModels(butler, patchIdList, bgModelList)
            pool.mapToPrevious(self.addNextVisit, extractPatchData(visit), bgModel)

        pool.mapToPrevious(self.finalize, patchRefList)

    def firstVisit(self, cache, patchData, visit):
        """Initialise background reference with the first visit

        This method is intended to be run by slave nodes under the process pool.

        Because we can only iterate over a single object, the
        'patchData' includes multiple elements:
            * patchRef: data reference for patch
            * calexpRefList: List of calexp data references overlapping patch

        @param cache: Process pool cache
        @param patchData: Input data for patch
        @param visit: Visit name
        """
        patchRef = patchData.patchRef
        calexpRefList = patchData.calexpRefList
        patchIndex = getPatchIndex(patchRef.dataId)
        skyInfo = getSkyInfo(cache.coaddName, patchRef)
        cache.skyInfo = skyInfo
        dataId = dict(zip(cache.visitKeys, visit))
        if patchRef.datasetExists(cache.warpType, **dataId):
            self.log.info("%s: Creating background reference for %s from %s" %
                          (NODE, patchIndex, dataId))
            warp = self.getWarp(patchRef, dataId, cache.warpType)
            weight = self.generateWeight(skyInfo, calexpRefList)
            self.zeroOutBad(warp, weight)
            warp.getMaskedImage().__imul__(weight)
        else:
            # Create a blank image we can use as a base
            self.log.info("%s: Creating blank background reference for %s" % (NODE, patchIndex))
            image = afwImage.MaskedImageF(skyInfo.bbox)
            image.setXY0(skyInfo.bbox.getMin())
            maskVal = image.getMask().getPlaneBitMask(self.config.mask)
            image.set((0.0, maskVal, 0.0))
            warp = afwImage.makeExposure(image, skyInfo.wcs)
            weight = afwImage.ImageF(skyInfo.bbox)
            weight.set(0.0)
        cache.bgRef = warp
        cache.bgWeight = weight

    def getWarp(self, patchRef, dataId, datasetType="deepCoadd_tempExp"):
        warp = patchRef.get(datasetType, immediate=True, **dataId)
        self.scaling.computeImageScaler(warp, patchRef).scaleMaskedImage(warp.getMaskedImage())
        return warp

    def matchNextVisit(self, cache, patchRef, visit):
        """Match a visit to the existing background reference

        This method is intended to be run by slave nodes under the process pool.

        @param cache: Process pool cache
        @param patchRef: Input data for patch
        @param visit: Visit name
        @return background model
        """
        cache.visit = visit
        patchIndex = getPatchIndex(patchRef.dataId)
        dataId = dict(zip(cache.visitKeys, visit))
        self.log.info("%s: Matching %s to %s" % (NODE, patchIndex, dataId))
        if not patchRef.datasetExists(cache.warpType, **dataId):
            # Nothing to do
            cache.warp = None
            return None
        warp = self.getWarp(patchRef, dataId, cache.warpType)
        cache.warp = warp
        bgRef = cache.bgRef.clone()
        bgRef.getMaskedImage().__idiv__(cache.bgWeight)

        bgModel = self.matching.run(bgRef, warp, inPlace=True)
        self.zeroOutBad(bgModel.getStatsImage().getImage())

        return bgModel

    def zeroOutBad(self, image, *others):
        """Zero out bad and non-finite pixels in an image

        The same pixels are also zeroed out in the 'others'.
        """
        others = list(others)
        if hasattr(image, "getMaskedImage"):
            image = image.getMaskedImage()
        if hasattr(image, "getMask"):
            # A MaskedImage
            mask = image.getMask()
            bitmask = mask.getPlaneBitMask(self.config.mask)
            bad = numpy.logical_or(mask.getArray() & bitmask > 0,
                                   numpy.logical_not(numpy.isfinite(image.getImage().getArray())))
            others.append(image.getImage())
            others.append(image.getVariance())
        else:
            bad = numpy.logical_not(numpy.isfinite(image.getArray()))
            others.append(image)

        for im in others:
            array = im.getArray()
            array[bad] = 0.0

    def mergeBackgroundModels(self, butler, patchIdList, bgModelList, coaddName="deep"):
        """Merge background models from different patches

        This ensures all patches use a common background model for
        matching to the reference.

        This implementation is extremely inefficient, and could be
        improved by adding support for this operation in the
        BackgroundMI class.

        @param butler: Data butler
        @param patchIdList: List of patch data identifiers
        @param bgModelList: Corresponding list of patch background models
        @param coaddName: Name of coadd (for retrieving skymap)
        @return merged background model
        """
        skymap = None
        statsImage = None
        # XXX configure stats
        statCtrl = afwMath.StatisticsControl()
        statCtrl.setAndMask(afwImage.MaskU.getPlaneBitMask(self.config.mask))
        stat = afwMath.MEAN
        for patchId, bgModel in zip(patchIdList, bgModelList):
            if bgModel is None:
                continue
            if not skymap:
                # Set up the patch
                skymap = butler.get(coaddName + "Coadd_skyMap", patchId)
                tract = skymap[patchId["tract"]]

                tractBox = tract.getBBox()
                x0, y0 = tractBox.getMin()
                # Size of tract
                xSize, ySize = tractBox.getMax() - afwGeom.Extent2I(x0, y0)
                # Number of samples for background model
                xNum = int(math.ceil(float(xSize)/self.config.bgSize))
                yNum = int(math.ceil(float(ySize)/self.config.bgSize))
                xBounds = numpy.linspace(x0, xSize, xNum + 1, True).astype(int)
                yBounds = numpy.linspace(x0, ySize, yNum + 1, True).astype(int)

                statsNum = afwImage.ImageF(xNum, yNum) # float so we can easily divide later
                statsNum.set(0.0)
                statsImage = afwImage.ImageF(xNum, yNum)
                statsImage.set(0.0)

                def statsIter():
                    """Return coordinates on stats image, box on tract"""
                    for ix, iy in itertools.product(range(xNum), range(yNum)):
                        yield ix, iy, afwGeom.Box2I(afwGeom.Point2I(int(xBounds[ix]), int(yBounds[iy])),
                                                    afwGeom.Point2I(int(xBounds[ix+1]), int(yBounds[iy+1])))

            # XXX this is terribly inefficient
            bgImage = bgModel.getImageF("AKIMA_SPLINE", "REDUCE_INTERP_ORDER")
            patchIndex = getPatchIndex(patchId)
            patchBox = tract[patchIndex].getOuterBBox()
            bgImage.setXY0(patchBox.getMin())
            for ix, iy, box in statsIter():
                try:
                    # XXX shrink box if only just a bit big
                    subImage = afwImage.ImageF(bgImage, box, afwImage.PARENT)
                    value = afwMath.makeStatistics(subImage, stat, statCtrl).getValue()
                except:
                    continue
                statsImage.set(ix, iy, statsImage.get(ix, iy) + value)
                statsNum.set(ix, iy, statsNum.get(ix, iy) + 1)
            del bgImage

        if statsImage is None:
            raise RuntimeError("Unable to merge background models")

        statsImage /= statsNum
        statsImage = afwImage.makeMaskedImage(statsImage)
        statsMask = statsImage.getMask().getArray()
        bad = numpy.where(numpy.isnan(statsImage.getImage().getArray()))
        statsMask[bad] = statsImage.getMask().getPlaneBitMask("BAD")
        mergedModel = afwMath.BackgroundMI(tractBox, statsImage)
        self.zeroOutBad(mergedModel.getStatsImage().getImage())

        return mergedModel

    def addNextVisit(self, cache, patchData, bgModel):
        """Add visit to background reference

        This method is intended to be run by slave nodes under the process pool.

        Because we can only iterate over a single object, the
        'patchData' includes multiple elements:
            * patchRef: data reference for patch
            * calexpRefList: List of calexp data references overlapping patch

        This implementation is extremely inefficient, and could be
        improved by adding support for this operation in the
        BackgroundMI class.

        @param cache: Process pool cache
        @param patchData: Input data for patch
        @param bgModel: Merged tract background model
        """
        patchRef = patchData.patchRef
        calexpRefList = patchData.calexpRefList
        patchIndex = getPatchIndex(patchRef.dataId)
        self.log.info("%s: Adding visit %s to patch %s" % (NODE, cache.visit, patchIndex))
        skyInfo = cache.skyInfo
        # XXX this is awfully inefficient
        bgImage = bgModel.getImageF("AKIMA_SPLINE", "REDUCE_INTERP_ORDER")
        bgSubImage = afwImage.ImageF(bgImage, skyInfo.bbox)

        if self.debug:
            bgSubImage.writeFits("model-%s-%s.fits" % ("-".join(map(str, cache.visit)), "%d,%d" % patchIndex))

        if cache.warp is not None:
            warp = cache.warp
            warpImage = warp.getMaskedImage()
            warpImage += bgSubImage
            weight = self.generateWeight(skyInfo, calexpRefList)
            self.zeroOutBad(warp, weight)
        else:
            warpImage = bgSubImage
            if calexpRefList:
                weight = self.generateWeight(skyInfo, calexpRefList)
                self.zeroOutBad(warpImage, weight)
            else:
                weight = afwImage.ImageF(warpImage.getDimensions())
                weight.set(0.0)

        if self.debug:
            warpImage.writeFits("add-%s-%s.fits" % ("-".join(map(str, cache.visit)), "%d,%d" % patchIndex))

        warpImage *= weight
        cache.bgRef.getMaskedImage().__iadd__(warpImage)
        cache.bgWeight += weight

        # Set bad pixels in the background reference
        bgImage = cache.bgRef.getMaskedImage()
        mask = bgImage.getMask()
        maskArray = mask.getArray()
        bad = numpy.where(numpy.isnan(bgImage.getImage().getArray()))
        bitmask = mask.getPlaneBitMask(self.config.mask)
        maskArray &= ~bitmask
        maskArray[bad] = mask.getPlaneBitMask("BAD")

        if self.debug:
            bgRef = cache.bgRef.getMaskedImage().clone()
            bgRef /= cache.bgWeight
            bgRef.writeFits("bg-%s-%s.fits" % ("-".join(map(str, cache.visit)), "%d,%d" % patchIndex))

    def finalize(self, cache, patchRef):
        """Finish up background reference and write

        This method is intended to be run by slave nodes under the process pool.

        @param cache: Process pool cache
        @param patchRef: Patch data reference
        """
        bgRef = cache.bgRef
        bgRef.getMaskedImage().__idiv__(cache.bgWeight)
        patchRef.put(bgRef, cache.bgRefType)

    def generateWeight(self, skyInfo, calexpRefList):
        """Construct a weight

        Returns a lightweight version of the actual weight image;
        call 'getImage' for the actual image.

        @param skyInfo: struct with skymap information
        @param calexpRefList: list of calexp data references
        @return weight image
        """
        transformer = CoordinateTransformer(skyInfo.wcs, calexpRefList)
        boresight = transformer.getBoresight()
        weight = self.config.taper.apply(skyInfo.bbox, transformer.skyToPatch(boresight))
        if hasattr(weight, "getImage"):
            weight = weight.getImage()
        return weight



class BackgroundReferenceConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd of interest")
    visitKeys = ListField(dtype=str, default=["visit"], doc="dataId keys that identify an exposure")
    select = ConfigurableField(target=WcsSelectImagesTask, doc="Task to select input images")
    assign = ConfigurableField(target=AssignTask, doc="Task to assign inputs to patches")
    construct = ConfigurableField(target=ConstructionTask, doc="Task to construct background reference")

class BackgroundReferenceTask(PbsPoolTask):
    _DefaultName = "bgRef"
    ConfigClass = BackgroundReferenceConfig
    RunnerClass = CoaddTaskRunner

    def __init__(self, *args, **kwargs):
        super(BackgroundReferenceTask, self).__init__(*args, **kwargs)
        self.makeSubtask("select")
        self.makeSubtask("assign")
        self.makeSubtask("construct")

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
        assignData = self.assign.run(tract, dataRefList)

        self.construct.run(patchRefList, assignData.visits, assignData.overlaps,
                           coaddName=self.config.coaddName, visitKeys=self.config.visitKeys)

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

