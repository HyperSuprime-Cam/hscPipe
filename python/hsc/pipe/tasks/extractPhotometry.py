import numpy

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

from lsst.pex.config import Config, Field, DictField, ConfigurableField, ChoiceField
from lsst.pipe.base import Task, CmdLineTask, ArgumentParser, TaskRunner, DatasetArgument, Struct
from lsst.pipe.tasks.coaddBase import ExistingCoaddDataIdContainer

"""
Module for extracting photometry from forcedPhotCoadd or processCoadd results.

By default, will extract PSF mags from forcedPhotCoadd.

For processCoadd, add the overrides:
    root.srcPostfix = "src"
    from hsc.pipe.tasks.extractPhotometry import MatchJoinTask
    root.join.retarget(MatchJoinTask)
"""

class SubsetConfig(Config):
    onlyStars = Field(dtype=bool, default=True, doc="Select only stars?")
    magLimit = Field(dtype=float, optional=True, doc="Magnitude limit")
    band = Field(dtype=str, optional=True, doc="Band for magnitude limit")

    def validate(self):
        super(SubsetConfig, self).validate()
        if (self.magLimit is not None and self.band is None) or \
                (self.magLimit is None and self.band is not None):
            raise RuntimeError("magLimit and band must both be set")

class SubsetTask(Task):
    ConfigClass = SubsetConfig

    def run(self, catalogDict, calibDict, refFilter):
        return catalogDict

class IndexJoiner(object):
    """A functor that will join (in a database sense) two arrays

    We would like to say:
        leftIndices, rightIndices = join(leftIdList, rightIdList)
        left[leftIndices] = right[rightIndices]
    but numpy returns a *copy* rather than a view when doing fancy indexing,
    which means it's not so simple.  Instead, we provide a functor that will
    act on 'right' to produce an array the same length as 'left'.  So the
    above becomes:
        join = IndexJoiner(leftIdList, rightIdList)
        left[:] = join(right)
    """

    def __init__(self, left, right):
        """Constructor

        Indices must be sorted ascending and elements unique

        @param left: Indices for left side
        @param right: Indices for right side
        """
        numLeft = len(left)
        numRight = len(right)
        select = numpy.zeros(numLeft, dtype=bool)
        indices = numpy.zeros(numLeft, dtype=int)
        iRight = 0
        for iLeft, elemLeft in enumerate(left):
            while iRight < numRight and right[iRight] < elemLeft:
                iRight += 1
            if iRight == numRight:
                break
            if right[iRight] == elemLeft:
                select[iLeft] = True
                indices[iLeft] = iRight

        self.num = numLeft
        self.select = select
        self.indices = indices

    def __call__(self, right, blank=numpy.nan):
        """Produce an array from 'right' the same length as 'left'

        Missing elements are set to 'blank'
        """
        null = numpy.empty(self.num)
        null.fill(blank)
        return numpy.where(self.select, right[self.indices], null)

class CoordJoiner(IndexJoiner):
    """Joiner that works by matching coordinates"""
    def __init__(self, left, right, radius):
        """Constructor

        @param left: Coord array for left side
        @param right: Coord array for right side
        @param radius: Angle specifying maximum match radius
        """
        num = len(left)
        schema = afwTable.SimpleTable.makeMinimalSchema()
        idKey = schema.find("id").key
        coordKey = schema.find("coord")

        def create(coordList):
            cat = afwTable.SimpleCatalog(schema)
            cat.reserve(len(coordList))
            for i, coord in enumerate(coordList):
                record = cat.addNew()
                record.set(idKey, i)
                record.setCoord(coord)
            return cat

        left = create(left)
        right = create(right)

        matches = sorted(afwTable.matchRaDec(left, right, radius),
                         cmp=lambda x,y: cmp(x.first.get(idKey), y.first.get(idKey)))
        left = numpy.array([m.first.get(idKey) for m in matches])
        right = numpy.array([m.second.get(idKey) for m in matches])
        numMatches = len(matches)
        select = numpy.zeros(num, dtype=bool)
        indices = numpy.zeros(num, dtype=int)
        iMatch = 0
        for iLeft in range(num):
            while iMatch < numMatches and left[iMatch] < iLeft:
                iMatch += 1
            if iMatch == numMatches:
                break
            if left[iMatch] == iLeft:
                select[iLeft] = True
                indices[iLeft] = right[iMatch]

        self.num = num
        self.select = select
        self.indices = indices


class JoinConfig(Config):
    fluxName = Field(dtype=str, default="flux.psf", doc="Flux column name to extract")
    extendedName = Field(dtype=str, default="classification.extendedness", doc="Extended column name")

class JoinTask(Task):
    ConfigClass = JoinConfig
    def run(self, catalogDict, calibDict, refFilter):
        catalog, keys = self.createCatalog(catalogDict[refFilter], catalogDict.keys())

        for filterName, cat in catalogDict.items():
            join = self.getJoin(catalogDict[refFilter], cat)

            flux = cat[self.config.fluxName]
            fluxErr = cat[self.config.fluxName + ".err"]
            calib = calibDict[filterName]

            mag = numpy.array([calib.getMagnitude(f, e) for f,e in zip(flux, fluxErr)])
            catalog[getattr(keys, filterName)][:] = join(mag[:,0])
            catalog[getattr(keys, filterName + "Err")][:] = join(mag[:,1])

        return catalog

    def createCatalog(self, refCatalog, filterList):
        schema = afwTable.Schema()
        idKey = schema.addField("id", type="L", doc="Object identifier")
        raKey = schema.addField("ra", type=float, doc="Right Ascension (ICRS; degrees)")
        decKey = schema.addField("dec", type=float, doc="Declination (ICRS; degrees)")
        starGalKey = schema.addField("starNotGal", type="I", doc="True if believed to be stellar")
        magKeys = dict([(f, schema.addField(f, type=numpy.float32, doc="Magnitude in %s" % f))
                       for f in filterList] +
                       [(f + "Err", schema.addField(f + "_err", type=numpy.float32,
                                                    doc="Magnitude error in %s" % f))
                        for f in filterList])
        catalog = afwTable.BaseCatalog(schema)

        num = len(refCatalog)
        catalog.reserve(num)
        for i in range(num):
            catalog.addNew()

        for key in magKeys.values():
            catalog[key][:] = numpy.nan

        keys = Struct(id=idKey, ra=raKey, dec=decKey, starGal=starGalKey, **magKeys)
        self.copyReferences(catalog, keys, refCatalog)
        return catalog, keys

    def copyReferences(self, catalog, keys, refCatalog):
        catalog[keys.id][:] = refCatalog['id']
        catalog[keys.ra][:] = numpy.degrees(refCatalog['coord.ra'])
        catalog[keys.dec][:] = numpy.degrees(refCatalog['coord.dec'])
        catalog[keys.starGal][:] = numpy.logical_not(refCatalog[self.config.extendedName])

    def getJoin(self, refCatalog, otherCatalog):
        return IndexJoiner(refCatalog["objectId"], otherCatalog["objectId"])


class MatchJoinConfig(JoinConfig):
    radius = Field(dtype=float, default=0.5, doc="Matching radius (arcsec)")

class MatchJoinTask(JoinTask):
    ConfigClass = MatchJoinConfig

    def getJoin(self, refCatalog, otherCatalog):
        refCoords = [s.getCoord() for s in refCatalog]
        otherCoords = [s.getCoord() for s in otherCatalog]
        return CoordJoiner(refCoords, otherCoords, self.config.radius*afwGeom.arcseconds)


class ExtractPhotometryConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Coadd name with measuremements")
    srcPostfix = Field(dtype=str, default="forced_src", doc="Postfix for source dataset")
    filters = DictField(keytype=str, itemtype=str, default={}, doc="Filters to extract, and what to call them")
    subset = ConfigurableField(target=SubsetTask, doc="How to subset the catalog")
    join = ConfigurableField(target=JoinTask, doc="How to join the catalogues")

class ExtractPhotometryTask(CmdLineTask):
    ConfigClass = ExtractPhotometryConfig
    _DefaultName = "extractPhotometry"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=ExistingCoaddDataIdContainer)
        return parser

    def __init__(self, *args, **kwargs):
        super(ExtractPhotometryTask, self).__init__(*args, **kwargs)
        self.makeSubtask("subset")
        self.makeSubtask("join")

    def run(self, dataRef):
        if not "filter" in dataRef.dataId:
            raise RuntimeError("No filter specified in dataRef: %s" % (dataRef.dataId,))
        refFilter = self.config.filters[dataRef.dataId["filter"]]
        catalogDict, calibDict = self.read(dataRef)
        catalogDict = self.subset.run(catalogDict, calibDict, refFilter)
        catalog = self.join.run(catalogDict, calibDict, refFilter)
        self.write(catalog, dataRef)

    def read(self, dataRef):
        catalogDict = {}
        calibDict = {}
        for filterName, magName in self.config.filters.items():
            catalogDict[magName], calibDict[magName] = self.readCatalog(dataRef, filter=filterName)
        return catalogDict, calibDict

    def readCatalog(self, dataRef, **dataId):
        dataset = self.config.coaddName + "Coadd_" + self.config.srcPostfix
        catalog = dataRef.get(dataset, immediate=True, **dataId)
        calib = self.getCalib(dataRef, **dataId)
        return catalog, calib

    def getCalib(self, dataRef, **dataId):
        coaddDataset = self.config.coaddName + "Coadd"

        bbox = afwGeom.Box2I(afwGeom.Point2I(), afwGeom.Extent2I(1, 1))
        exp = dataRef.get(coaddDataset + "_sub", bbox=bbox, imageOrigin="LOCAL", immediate=True, **dataId)
        calib = exp.getCalib()
        calib.setThrowOnNegativeFlux(False)
        return calib

    def write(self, catalog, dataRef):
        dataRef.put(catalog, self.config.coaddName + "Coadd_extract")

    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass
