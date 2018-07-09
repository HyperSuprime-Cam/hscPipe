from __future__ import absolute_import, division, print_function

import os

from lsst.utils import getPackageDir
from lsst.afw.table import SOURCE_IO_NO_FOOTPRINTS
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.meas.base import ForcedPhotCoaddConfig, ForcedPhotCoaddTask
from lsst.pipe.base import ButlerInitializedTaskRunner
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner


__all__ = ["UndeblendedConvolvedConfig",
           "UndeblendedConvolvedTask",
           "UndeblendedConvolvedDriverConfig",
           "UndeblendedConvolvedDriverTask",
           ]


def doUndeblended(config, algName, fluxList=None):
    """Activate undeblended measurements of algorithm

    Parameters
    ----------
    algName : `str`
        Algorithm name.
    fluxList : `list` of `str`, or `None`
        List of flux columns to register for aperture correction. If `None`,
        then this will be the `algName` appended with `_flux`.
    """
    if algName not in config.measurement.plugins:
        raise RuntimeError("Can't find algorithm %s" % (algName,))
    if fluxList is None:
        fluxList = [algName + "_flux"]
    config.measurement.undeblended.names.add(algName)
    config.measurement.undeblended[algName] = config.measurement.plugins[algName]
    for flux in fluxList:
        config.applyApCorr.proxies["undeblended_" + flux] = flux


class UndeblendedConvolvedConfig(ForcedPhotCoaddConfig):
    """Configuration for forced convolved measurement of undeblended sources"""
    def setDefaults(self):
        """Set defaults for configuration

        Disable all measurements we don't care about, and activate the convolved
        measurements for both deblended and undeblended sources.
        """
        ForcedPhotCoaddConfig.setDefaults(self)

        self.measurement.slots.psfFlux = None
        self.measurement.slots.instFlux = None
        self.measurement.slots.shape = None
        self.measurement.slots.centroid = "base_TransformedCentroid"

        self.measurement.load(os.path.join(getPackageDir("obs_subaru"), "config", "convolvedFluxes.py"))
        doUndeblended(self, "ext_convolved_ConvolvedFlux",
                      self.measurement.plugins["ext_convolved_ConvolvedFlux"].getAllResultNames())
        self.measurement.undeblended["ext_convolved_ConvolvedFlux"].registerForApCorr = False
        # Disable all regular deblended convolved measurement except the essential centroid
        self.measurement.plugins.names = ["base_TransformedCentroid"]


class UndeblendedConvolvedTask(ForcedPhotCoaddTask):
    """Task for forced convolved measurement of undeblended sources"""
    ConfigClass = UndeblendedConvolvedConfig
    _DefaultName = "undeblendedConvolved"

    def writeOutput(self, dataRef, sources):
        """Write forced source table

        We don't want to write over the forced measurements written as part of
        the regular pipeline, so we need to write the results in a different
        file.

        Parameters
        ----------
        dataRef : `lsst.daf.persistence.ButlerDataRef`
            Butler data reference.
        sources : `lsst.afw.table.SourceCatalog`
            Catalog of sources to write.
        """
        dataRef.put(sources, self.dataPrefix + "forced_src_undeblendedConvolved",
                    flags=SOURCE_IO_NO_FOOTPRINTS)

    def _getConfigName(self):
        """Return name of configuration dataset"""
        return "undeblendedConvolved_config"

    def _getMetadataName(self):
        """Disable metadata writing"""
        return None


class UndeblendedConvolvedRunner(BatchTaskRunner, ButlerInitializedTaskRunner):
    """Run batches, and initialize Task using a butler"""
    pass


class UndeblendedConvolvedDriverConfig(Config):
    """Configuration for parallel forced convolved measurement of undeblended
    sources
    """
    undeblendedConvolved = ConfigurableField(target=UndeblendedConvolvedTask, doc="Main task")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd type")


class UndeblendedConvolvedDriverTask(BatchParallelTask):
    """Task for parallel forced convolved measurement of undeblended sources"""
    ConfigClass = UndeblendedConvolvedDriverConfig
    _DefaultName = "undeblendedConvolvedDriver"
    RunnerClass = UndeblendedConvolvedRunner

    def __init__(self, butler, *args, **kwargs):
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.makeSubtask("undeblendedConvolved", butler=butler)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        kwargs.pop("add_help", False)
        return UndeblendedConvolvedTask._makeArgumentParser(*args, **kwargs)

    def run(self, dataRef):
        with self.logOperation("processing %s" % (dataRef.dataId,)):
            return self.undeblendedConvolved.run(dataRef)

    def _getMetadataName(self):
        """Disable metadata writing"""
        return None
