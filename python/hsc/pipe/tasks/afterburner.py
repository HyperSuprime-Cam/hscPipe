#!/usr/bin/env python

"""
Afterburner for the 2016 May-June production run.

We will run the following:
* PSF-matched aperture fluxes on coadd, for both parent and child
* Detection using the over-subtracted background and flag sources not detected
"""

import itertools
import contextlib

from lsst.pex.config import Config, Field, ConfigurableField, ListField
from lsst.pipe.base import Task, CmdLineTask, ArgumentParser, ButlerInitializedTaskRunner, Struct
from lsst.meas.algorithms import SourceDetectionTask, MeasureSourcesBuilder
from lsst.meas.algorithms.replaceWithNoise import ReplaceWithNoiseTask
from lsst.meas.algorithms.algorithmRegistry import AlgorithmRegistry
from lsst.meas.algorithms.measurement import ExposureContext, SourceContext
from lsst.meas.algorithms.measureApCorr import KeyTuple, MeasureApCorrTask
from lsst.pipe.tasks.coaddBase import ExistingCoaddDataIdContainer
from hsc.pipe.base.parallel import BatchCmdLineTask, abortOnError
from hsc.pipe.base.pool import Pool, startPool

import lsst.afw.table as afwTable
import lsst.afw.image as afwImage

import lsst.meas.extensions.convolved  # register measurement algorithm

class JunkConfig(Config):
    detection = ConfigurableField(target=SourceDetectionTask, doc="Detect sources")
    matchRadius = Field(dtype=float, default=0.25, doc="Matching radius (arcsec)")

    def setDefaults(self):
        Config.setDefaults(self)
        # Same overrides as for detectCoaddSources, but we'll turn off background subtraction and turn on
        # the footprint background
        self.detection.thresholdType = "pixel_stdev"
        self.detection.isotropicGrow = True
        self.detection.reEstimateBackground = False  # Don't do a lasting background subtraction on the image
        self.detection.doFootprintBackground = True  # Do a temporary background subtraction for detection

class JunkTask(Task):
    """Identify junk reference sources

    We run detection on the exposure using the temporary background
    subtraction (junk suppression), match with the reference catalog,
    and flag sources that are detected.
    """
    ConfigClass = JunkConfig

    def __init__(self, schema, **kwargs):
        """Ctor

        @param schema: output catalog schema (updated with detection flags)
        """
        Task.__init__(self, **kwargs)
        self.makeSubtask("detection")  # Not providing schema; we don't care about the schema for detection
        self.detectedKey = schema.addField("detected.notJunk", type="Flag",
                                           doc="Detected after junk supression")

    def run(self, exposure, catalog):
        """Run the junk identification

        @param exposure: exposure on which to run detection
        @param catalog: reference catalog
        """
        exposure = exposure.Factory(exposure, True)  # copy exposure so it's not clobbered (esp. the mask)
        schema = afwTable.SourceTable.makeMinimalSchema()  # throwaway schema for detection
        centroidKey = schema.addField("centroid", type="PointD", doc="Centroid position")
        table = afwTable.SourceTable.make(schema)
        table.defineCentroid(centroidKey)
        results = self.detection.makeSourceCatalog(table, exposure)
        sources = results.sources
        children = afwTable.SourceCatalog(schema)
        numPeaks = sum(len(src.getFootprint().getPeaks()) for src in sources)
        self.log.info("Detected %d parents with %d children" % (len(sources), numPeaks))
        children.reserve(numPeaks)
        children.table.defineCentroid(centroidKey)
        for src in sources:
            peaks = src.getFootprint().getPeaks()
            for pk in peaks:
                child = children.addNew()
                child.set(centroidKey, pk.getCentroid())
        radius = self.config.matchRadius/exposure.getWcs().pixelScale().asArcseconds()
        matches = afwTable.matchXy(catalog, children, radius)
        self.log.info("Matched %d children compared with reference %d objects" %
                      (len(matches), len(catalog)))
        for mm in matches:
            mm.first.set(self.detectedKey, True)

@contextlib.contextmanager
def dummyContext(*args, **kwargs):
    """A do-nothing context manager

    Useful for NOT replacing sources with noise.
    """
    try:
        yield
    except Exception as exc:  # Don't catch KeyboardInterrupt
        print "Swallowed exception: %s" % (exc,)
        pass


class AfterburnerApCorrTask(MeasureApCorrTask):
    def __init__(self, schema, columns, **kwargs):
        """Ctor

        Different setup than for MeasureApCorrTask because we get the list of columns to correct
        from a config parameter rather than a registry. The registry won't work for our purposes
        because we want to do multiple corrections under flux.convolved.
        """
        lsst.pipe.base.Task.__init__(self, **kwargs)
        self.reference = KeyTuple(self.config.reference, schema)
        self.toCorrect = {col: KeyTuple(col, schema) for col in columns}
        self.inputFilterFlag = schema.find(self.config.inputFilterFlag).key
        self.keys = {col: schema.addField(col + ".apCorr", type="D", doc="Aperture correction for " + col) for
                     col in columns}
        
    def run(self, bbox, catalog):
        """Measure and apply aperture correction"""
        apCorrMap = MeasureApCorrTask.run(self, bbox, catalog)
        self.apply(catalog, apCorrMap)
        return apCorrMap

    def apply(self, catalog, apCorrMap):
        """Apply aperture correction"""
        xx, yy = catalog["centroid.sdss.x"], catalog["centroid.sdss.y"]
        for col in self.toCorrect:
            corr = apCorrMap[col].evaluate(xx, yy)
            catalog[self.keys[col]][:] = corr
            catalog[col][:] *= corr


class MeasureConfig(Config):
    algorithms = AlgorithmRegistry.all.makeField(multi=True, default=["flux.naive", "flux.convolved"],
                                                 doc="Measurement algorithms that will be run")
    replaceWithNoise = ConfigurableField(target=ReplaceWithNoiseTask,
                                         doc="Replacing other sources by noise when measuring sources")
    measureApCorr = ConfigurableField(target=AfterburnerApCorrTask, doc="Aperture correction")
    apCorrColumns = ListField(dtype=str, doc="Columns for aperture correction",
                              default=(list("flux.convolved.%d.%d" % xx for
                                            xx in itertools.product(range(3), range(3))) +
                                       list("flux.convolved.%d.kron" % xx for xx in range(3))))
    doReplaceWithNoise = True  # Will be used internally

    def setDefaults(self):
        Config.setDefaults(self)
        self.algorithms["flux.naive"].radius = 12.0

        self.measureApCorr.reference = "child.flux.naive"
        # More extreme clipping because we expect the distribution to have longer tails on the coadd
        self.measureApCorr.numSigmaClip = 2.5
        self.measureApCorr.numIter = 10

class MeasureFamiliesTask(Task):
    """Like the traditional MeasurementTask, but measures parents and children"""
    ConfigClass = MeasureConfig

    def __init__(self, schema, **kwargs):
        Task.__init__(self, **kwargs)
        self.makeSubtask("replaceWithNoise")
        self.childMeasurer = self.makeMeasureSources("child.", schema, self.metadata)
        self.parentMeasurer = self.makeMeasureSources("parent.", schema, self.metadata)
        apCorrColumns = ["%s.%s" % xx for
                         xx in itertools.product(("child", "parent"), self.config.apCorrColumns)]
        self.makeSubtask("measureApCorr", schema=schema, columns=apCorrColumns)

    def makeMeasureSources(self, prefix, schema, metadata):
        builder = MeasureSourcesBuilder(prefix, True)
        builder.addAlgorithms(self.config.algorithms.apply())
        return builder.build(schema, metadata)

    def run(self, exposure, sources, references):
        """Measure on children and parents

        @param exposure: exposure on which to measure
        @param sources: output source catalog; will be modified with measurements
        @param references: reference catalog
        """
        bbox = exposure.getBBox(afwImage.PARENT)
        self.log.info("Measuring on children")
        self.measure(self.childMeasurer, ExposureContext, SourceContext,
                     exposure, sources, references)
        self.log.info("Measuring on parents")
        self.measure(self.parentMeasurer, dummyContext, dummyContext,  # no noise replacement
                     exposure, sources, references)
        apCorrMap = self.measureApCorr.run(bbox, sources)
        return Struct(apCorrMap=apCorrMap)

    # Do-nothing methods
    #
    # These are required by the ReplaceWithNoiseTask. In the original SourceMeasurementTask, these are used
    # only for display/debugging (i.e., have no impact on the measurements themselves).
    def preMeasureHook(*args, **kwargs): pass
    def postMeasureHook(*args, **kwargs): pass
    def preSingleMeasureHook(*args, **kwargs): pass
    def postSingleMeasureHook(*args, **kwargs): pass

    def measure(self, measurer, expContext, srcContext, exposure, sources, references):
        """Measure sources

        @param measurer: MeasureSources object
        @param expContext: exposure context; used to replace all objects with noise
        @param srcContext: source context; used to restore object of interest for measurement
        @param exposure: exposure on which to measure
        @param sources: output source catalog; will be modified with measurements
        @param references: reference catalog
        """
        with expContext(self, exposure, sources):  # Replace all images with noise
            for i, (source, ref) in enumerate(zip(sources, references)):
                with srcContext(self, exposure, sources, i):  # Restore source of interest
                    measurer.applyForced(source, exposure, ref, exposure.getWcs())


class AfterburnerConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    junk = ConfigurableField(target=JunkTask, doc="Identify junk")
    measureFamilies = ConfigurableField(target=MeasureFamiliesTask, doc="Measure families")


class AfterburnerTask(CmdLineTask):
    ConfigClass = AfterburnerConfig
    _DefaultName = "afterburner"
    RunnerClass = ButlerInitializedTaskRunner

    def __init__(self, butler, **kwargs):
        CmdLineTask.__init__(self, **kwargs)
        schema = butler.get(self.config.coaddName + "Coadd_ref_schema", immediate=True).schema
        self.mapper = afwTable.SchemaMapper(schema)
        self.mapper.addMinimalSchema(afwTable.SourceTable.makeMinimalSchema())
        self.mapper.addMapping(schema["centroid.sdss"].asKey())
        self.mapper.addMapping(schema["calib.psf.used"].asKey())
        self.mapper.addMapping(schema["calib.psf.candidate"].asKey())
        self.schema = self.mapper.getOutputSchema()
        self.makeSubtask("junk", schema=self.schema)
        self.makeSubtask("measureFamilies", schema=self.schema)

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd_calexp",
                               ContainerClass=ExistingCoaddDataIdContainer,
                               help="data ID, e.g. --id tract=12345 patch=1,2 filter=g^r^i")
        return parser

    def run(self, dataRef):
        references = self.readReferences(dataRef)
        measurements = dataRef.get(self.config.coaddName + "Coadd_meas", immediate=True)
        self.log.info("Read %d reference sources" % (len(references),))

        sources = self.makeCatalog(references, measurements)
        exposure = dataRef.get(self.config.coaddName + "Coadd_calexp", immediate=True)

        self.junk.run(exposure, sources)
        results = self.measureFamilies.run(exposure, sources, references)
        dataRef.put(sources, self.config.coaddName + "Coadd_afterburner",
                    flags=afwTable.SOURCE_IO_NO_FOOTPRINTS)

        filename = dataRef.get(self.config.coaddName + "Coadd_afterburnerApCorr_filename")[0]
        results.apCorrMap.writeFits(filename)

    def readReferences(self, dataRef):
        references = dataRef.get(self.config.coaddName + "Coadd_ref", immediate=True)
        if references.schema == self.mapper.getInputSchema():
            return references
        # The May/June production run added columns with i2 that aren't in the Dec production run
        schema = self.mapper.getInputSchema()
        mapper = afwTable.SchemaMapper(references.schema)
        for ff in schema:
            name = ff.field.getName()
            if name in references.schema:
                mapper.addMapping(references.schema[name].asKey(), schema[name].asField())
            else:
                mapper.addOutputField(ff.field)

        replacement = afwTable.SourceCatalog(schema)
        replacement.extend(references, mapper)
        replacement.table.defineCentroid(references.table.getCentroidDefinition())
        replacement.table.defineShape(references.table.getShapeDefinition())
        replacement.table.definePsfFlux(references.table.getPsfFluxDefinition())
        replacement.table.defineModelFlux(references.table.getModelFluxDefinition())
        replacement.table.defineApFlux(references.table.getApFluxDefinition())
        replacement.table.defineInstFlux(references.table.getInstFluxDefinition())
        replacement.table.defineCalibFlux(references.table.getCalibFluxDefinition())
        return replacement

    def makeCatalog(self, references, measurements):
        sources = afwTable.SourceCatalog(self.schema)
        sources.table.defineCentroid("centroid.sdss")
        sources.extend(references, mapper=self.mapper)
        for src, meas in zip(sources, measurements):
            assert src.getId() == meas.getId()
            src.setFootprint(meas.getFootprint())

        return sources

    # XXX for testing only
    def _getConfigName(self):
        return None
    def _getMetadataName(self):
        return None
    def _getEupsVersionsName(self):
        return None


class AfterburnerBatchTaskRunner(ButlerInitializedTaskRunner):
    """Run a Task individually on a list of inputs using the MPI process pool"""
    def __init__(self, *args, **kwargs):
        """Constructor

        Warn if the user specified multiprocessing.
        """
        ButlerInitializedTaskRunner.__init__(self, *args, **kwargs)
        if self.numProcesses > 1:
            self.log.warn("Multiprocessing arguments (-j %d) ignored since using batch processing" %
                          self.numProcesses)
            self.numProcesses = 1

    def run(self, parsedCmd):
        """Run the task on all targets

        Sole input is the result of parsing the command-line with the ArgumentParser.

        Output is None if 'precall' failed; otherwise it is a list of calling ourself
        on each element of the target list from the 'getTargetList' method.
        """
        resultList = None

        self.prepareForMultiProcessing()
        pool = Pool()

        if self.precall(parsedCmd):
            targetList = self.getTargetList(parsedCmd)
            if len(targetList) > 0:
                parsedCmd.log.info("Processing %d targets with a pool of %d processes..." %
                                   (len(targetList), pool.size))
                # Run the task using self.__call__
                resultList = pool.map(self, targetList)
            else:
                parsedCmd.log.warn("Not running the task because there is no data to process; "
                    "you may preview data using \"--show data\"")
                resultList = []

        return resultList

    @abortOnError
    def __call__(self, cache, args):
        """Run the Task on a single target

        Strips out the process pool 'cache' argument.

        'args' are those arguments provided by the getTargetList method.

        Brings down the entire job if an exception is not caught (i.e., --doraise).
        """
        return ButlerInitializedTaskRunner.__call__(self, args)

class AfterburnerBatchTask(BatchCmdLineTask, AfterburnerTask):
    RunnerClass = AfterburnerBatchTaskRunner

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        """Build an ArgumentParser

        Removes the batch-specific parts in order to delegate to the parent classes.
        """
        kwargs.pop("doBatch", False)
        kwargs.pop("add_help", False)
        return super(BatchCmdLineTask, cls)._makeArgumentParser(*args, **kwargs)

    @classmethod
    def parseAndRun(cls, *args, **kwargs):
        """Parse an argument list and run the command

        This is the entry point when we run in earnest, so start the process pool
        so that the worker nodes don't go any further.
        """
        pool = startPool()
        results = super(AfterburnerBatchTask, cls).parseAndRun(*args, **kwargs)
        pool.exit()
        return results
