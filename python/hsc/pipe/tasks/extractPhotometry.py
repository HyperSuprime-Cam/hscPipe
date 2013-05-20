import numpy

import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

from lsst.pex.config import Config, Field, DictField, ConfigurableField, FieldValidationError
from lsst.pipe.base import Task, CmdLineTask, ArgumentParser, TaskRunner
from lsst.pipe.tasks.coaddBase import ExistingCoaddDataIdContainer

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

    def run(self, catalog):
        return catalog


class ExtractPhotometryConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Coadd name with measurements")
    filters = DictField(keytype=str, itemtype=str, default={}, doc="Filters to extract, and what to call them")
    fluxName = Field(dtype=str, default="flux.psf", doc="Flux column name to extract")
    extendedName = Field(dtype=str, default="classification.extendedness", doc="Extended column name")
    subset = ConfigurableField(target=SubsetTask, doc="How to subset the catalog")

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

    def run(self, dataRef):
        catalog = self.read(dataRef)
        catalog = self.subset.run(catalog)
        self.write(catalog, dataRef)

    def read(self, dataRef):
        coaddDataset = self.config.coaddName + "Coadd"

        # Reference list for forced photometry excludes sources in the boundary, so the length of the _src
        # and _forced_src catalogues doesn't match; for now, do forced photometry on the reference as well
        # to make the catalogues match.
        refDataset = coaddDataset + "_forced_src"
        forcedDataset = coaddDataset + "_forced_src"

        catalogDict = {}
        calibDict = {}

        refFilter = self.config.filters[dataRef.dataId["filter"]]
        catalogDict[refFilter] = dataRef.get(refDataset, immediate=True)
        calibDict[refFilter] = self.getCalib(dataRef, coaddDataset)

        for filterName, magName in self.config.filters.items():
            catalogDict[magName] = dataRef.get(forcedDataset, filter=filterName)
            calibDict[magName] = self.getCalib(dataRef, coaddDataset, filter=filterName)
            assert len(catalogDict[magName]) == len(catalogDict[refFilter]), "Catalog length mismatch"

        return self.join(catalogDict, calibDict, refFilter)

    def getCalib(self, dataRef, coaddDataset, **dataId):
        bbox = afwGeom.Box2I(afwGeom.Point2I(), afwGeom.Extent2I(1, 1))
        coadd = dataRef.get(coaddDataset + "_sub", bbox=bbox, imageOrigin="LOCAL", immediate=True, **dataId)
        calib = coadd.getCalib()
        calib.setThrowOnNegativeFlux(False)
        return calib

    def join(self, catalogDict, calibDict, refFilter):
        schema = afwTable.Schema()
        idKey = schema.addField("id", type="L", doc="Object identifier")
        raKey = schema.addField("ra", type=float, doc="Right Ascension (ICRS; degrees)")
        decKey = schema.addField("dec", type=float, doc="Declination (ICRS; degrees)")
        starGalKey = schema.addField("starNotGal", type="I", doc="True if believed to be stellar")
        magKeys = dict((f, schema.addField(f, type=numpy.float32, doc="Magnitude in %s" % f))
                       for f in catalogDict)
        errKeys = dict((f, schema.addField(f + "_err", type=numpy.float32, doc="Magnitude error in %s" % f))
                       for f in catalogDict)

        num = len(catalogDict[refFilter])
        catalog = afwTable.BaseCatalog(schema)
        catalog.reserve(num)
        for i in range(num):
            catalog.addNew()

        catalog[idKey][:] = catalogDict[refFilter]['id']
        catalog[raKey][:] = numpy.degrees(catalogDict[refFilter]['coord.ra'])
        catalog[decKey][:] = numpy.degrees(catalogDict[refFilter]['coord.dec'])
        catalog[starGalKey][:] = numpy.logical_not(catalogDict[refFilter][self.config.extendedName])

        for filterName, cat in catalogDict.items():
            flux = cat[self.config.fluxName]
            fluxErr = cat[self.config.fluxName + ".err"]
            calib = calibDict[filterName]

            magList = [calib.getMagnitude(f, e) for f,e in zip(flux, fluxErr)]
            catalog[magKeys[filterName]][:] = numpy.array([mag[0] for mag in magList])
            catalog[errKeys[filterName]][:] = numpy.array([mag[1] for mag in magList])

        return catalog

    def write(self, catalog, dataRef):
        dataRef.put(catalog, self.config.coaddName + "Coadd_extract")

    def writeConfig(self, *args, **kwargs):
        pass
    def writeSchema(self, *args, **kwargs):
        pass
    def writeMetadata(self, dataRef):
        pass
