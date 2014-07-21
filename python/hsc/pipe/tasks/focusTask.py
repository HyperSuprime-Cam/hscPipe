import os
import math
import collections

from lsst.pex.config import Config, Field, DictField, ConfigField, ConfigurableField
import lsst.afw.table as afwTable
import lsst.afw.cameraGeom as afwCG
import lsst.afw.cameraGeom.utils as afwCGU
import lsst.meas.algorithms as measAlg
from lsst.pipe.base import Task, Struct, ArgumentParser
from lsst.pipe.tasks.calibrate import CalibrateTask, InitialPsfConfig

from lsst.obs.subaru.isr import SubaruIsrTask

from hsc.pipe.base.pool import Pool
from hsc.pipe.base.parallel import BatchPoolTask
from hsc.pipe.base.butler import getDataRef

from .focus import FocusConfig, getDistanceFromFocus, haveSimpleShape

class ProcessFocusConfig(Config):
    focus = ConfigField(dtype=FocusConfig, doc="Focus determination")
    zemax = DictField(keytype=str, itemtype=str, default={},
                      doc="Mapping from filter name to zemax configuration filename")
    isr = ConfigurableField(target=SubaruIsrTask, doc="Instrument Signature Removal")
    background = ConfigField(dtype=measAlg.estimateBackground.ConfigClass, doc="Background removal")
    initialPsf = ConfigField(dtype=InitialPsfConfig, doc=InitialPsfConfig.__doc__)
    detection = ConfigurableField(target=measAlg.SourceDetectionTask, doc="Source detection")
    measurement = ConfigurableField(target=measAlg.SourceMeasurementTask, doc="Source measurement")
    starSelector = measAlg.starSelectorRegistry.makeField("Star selection algorithm", default="secondMoment")
    doWrite = Field(dtype=bool, default=True, doc="Write processed image?")

    def setDefaults(self):
        """These defaults are suitable for HSC, but may be useful
        for other cameras if the focus code is employed elsewhere.
        """
        super(ProcessFocusConfig, self).setDefaults()
        zemaxBase = os.path.join(os.environ["OBS_SUBARU_DIR"], "hsc", "zemax_config%d_0.0.dat")
        self.zemax = dict([(f, zemaxBase % n) for f,n in [('g', 9), ('r', 1), ('i', 3), ('z', 5), ('y', 7)]])
        self.load(os.path.join(os.environ["OBS_SUBARU_DIR"], "config", "hsc", "isr.py"))
        self.initialPsf.fwhm = 1.5 # arcsec
        self.initialPsf.size = 21 # pixels
        self.detection.includeThresholdMultiplier = 3.0
        self.measurement.centroider.name = "centroid.gaussian"
        self.measurement.slots.centroid = "centroid.gaussian"
        # set up simple shape, if available (because focus calibrations are for that)
        # If it's not available, we'll crash later; but we don't want to crash here (brings everything down)!
        if haveSimpleShape:
            self.measurement.algorithms.names.add("shape.simple")
            self.measurement.algorithms["shape.simple"].sigma = 5.0 # pixels
            self.measurement.slots.shape = "shape.simple"
        # set up background estimate
        self.background.ignoredPixelMask = ['EDGE', 'NO_DATA', 'DETECTED', 'DETECTED_NEGATIVE', 'BAD']
        self.detection.background.algorithm='LINEAR'
        self.starSelector.name = "objectSize"
        self.starSelector["objectSize"].badFlags = ["flags.pixel.edge",
                                                    "flags.pixel.interpolated.center",
                                                    "flags.pixel.saturated.center",
                                                    "flags.pixel.bad",
                                                    ]
        self.starSelector["objectSize"].sourceFluxField = "flux.gaussian"
        self.starSelector["objectSize"].widthMax = 20.0
        self.starSelector["objectSize"].widthStdAllowed = 5.0


class ProcessFocusTask(BatchPoolTask):
    ConfigClass = ProcessFocusConfig
    _DefaultName = "processFocus"
    installInitialPsf = CalibrateTask.__dict__["installInitialPsf"]

    def __init__(self, *args, **kwargs):
        super(ProcessFocusTask, self).__init__(*args, **kwargs)
        self.makeSubtask("isr")
        self.schema = afwTable.SourceTable.makeMinimalSchema()
        self.makeSubtask("detection", schema=self.schema)
        self.makeSubtask("measurement", schema=self.schema)
        self.starSelector = self.config.starSelector.apply()
        self.candidateKey = self.schema.addField(
            "calib.psf.candidate", type="Flag",
            doc=("Flag set if the source was a candidate for PSF determination, "
                 "as determined by the '%s' star selector.") % self.config.starSelector.name
        )

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numNodes, numProcs):
        config = parsedCmd.config
        numCcds = len(config.focus.aboveList) + len(config.focus.belowList)
        numCycles = int(math.ceil(numCcds/float(numNodes*numProcs)))
        numExps = len(cls.RunnerClass.getTargetList(parsedCmd))
        return time*numExps*numCycles

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        doBatch = kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="processFocus", *args, **kwargs)
        parser.add_id_argument("--id", datasetType="raw", level="visit",
                               help="data ID, e.g. --id visit=12345")
        return parser

    def run(self, expRef):
        """Measure focus for exposure

        This method is the top-level for running the focus measurement
        as a stand-alone BatchPoolTask.

        Only the master node runs this method.
        """
        pool = Pool("processFocus")
        pool.cacheClear()
        pool.storeSet(butler=expRef.getButler())

        dataIdList = sorted([ccdRef.dataId for ccdRef in expRef.subItems("ccd") if
                             ccdRef.datasetExists("raw") and self.isFocus(ccdRef)])

        results = pool.map(self.processPool, dataIdList)

        camera = expRef.get("camera")
        plotFilename = expRef.get("focusPlot_filename")
        focus = self.measureFocus(results, camera, plotFilename)
        self.log.info("Focus result for %s: %s" % (expRef.dataId, focus))
        return focus

    def isFocus(self, dataRef):
        """Is the provided dataRef for a focus CCD?"""
        ccdId = dataRef.dataId["ccd"]
        return self.config.focus.isFocusCcd(ccdId)

    def processPool(self, cache, dataId):
        """Process focus CCD under pool

        This is a mediator for the 'process' method when running
        under the Pool.

        Only slave nodes run this method.

        @param cache: Pool cache
        @param dataId: Data identifier for CCD
        @return Processing results (from 'process' method)
        """
        try:
            return self.process(getDataRef(cache.butler, dataId))
        except Exception as e:
            self.log.warn("Failed to process %s (%s): %s" % (dataId, e.__class__.__name__, e))
            return None

    def process(self, dataRef):
        """Process focus CCD in preparation for focus measurement

        @param dataRef: Data reference for CCD
        @return Struct(sources: source measurements,
                       ccdId: CCD number,
                       filterName: name of filter,
                       dims: exposure dimensions
                       )
        """
        import lsstDebug
        display = lsstDebug.Info(__name__).display

        exp = self.isr.run(dataRef).exposure

        if display:
            import lsst.afw.display.ds9 as ds9
            ds9.mtv(exp, title="Post-ISR", frame=1)

        self.installInitialPsf(exp)
        bg, exp = measAlg.estimateBackground(exp, self.config.background, subtract=True)

        if display:
            ds9.mtv(exp, title="Post-background", frame=2)

        table = afwTable.SourceTable.make(self.schema, afwTable.IdFactory.makeSimple())
        detRet = self.detection.makeSourceCatalog(table, exp)
        sources = detRet.sources
        self.measurement.run(exp, sources)
        psfCandidateList = self.starSelector.selectStars(exp, sources)
        if psfCandidateList:
            for cand in psfCandidateList:
                source = cand.getSource()
                source.set(self.candidateKey, True)

        if display:
            ds9.mtv(exp, title="Post-measurement", frame=3)
            with ds9.Buffering():
                for s in sources:
                    ds9.dot("o", s.getX(), s.getY(), frame=3,
                            ctype=ds9.GREEN if s.get("calib.psf.candidate") else ds9.RED)
            import pdb;pdb.set_trace() # pause to allow inspection

        filterName = exp.getFilter().getName()

        if self.config.doWrite:
            dataRef.put(sources, "icSrc")
            dataRef.put(exp, "visitim")

        return Struct(sources=sources, ccdId=dataRef.dataId["ccd"], filterName=filterName,
                      dims=exp.getDimensions())

    def measureFocus(self, resultsList, camera, plotFilename=None):
        """Measure focus from combining individual CCDs

        @param resultsList: Results of processing individual CCDs
        @param camera: Camera object
        @param plotFilename: Name of file for plot
        @return tuple(corrected distance from focus,
                      error in corrected distance from focus,
                      uncorrected distance from focus,
                      error in uncorrected distance from focus,
                      )
        """
        resultsList = [res for res in resultsList if res is not None]
        sources = dict((res.ccdId, res.sources) for res in resultsList)

        ccdList = [afwCG.cast_Ccd(ccd) for raft in camera for ccd in afwCG.cast_Raft(raft)]
        ccds = dict((ccd.getId().getSerial(), ccd) for ccd in ccdList)

        dims = dict((res.ccdId, res.dims) for res in resultsList)
        filterSet = set([res.filterName for res in resultsList])
        assert len(filterSet) == 1
        filterName = filterSet.pop()
        zemax = self.config.zemax[filterName]

        return getDistanceFromFocus(sources, ccds, dims, zemax, self.config.focus, plotFilename=plotFilename)
