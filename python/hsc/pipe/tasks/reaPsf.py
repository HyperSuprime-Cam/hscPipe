from lsst.pipe.base import CmdLineTask
from lsst.pex.config import Config, ConfigurableField
from lsst.pipe.tasks.measurePsf import MeasurePsfTask


class ReaPsfConfig(Config):
    measurePsf = ConfigurableField(target=MeasurePsfTask, doc="How to measure a PSF")

    def setDefaults(self):
        self.measurePsf.starSelector.name = "objectSize"
        self.measurePsf.starSelector["objectSize"].sourceFluxField = "initial.flux.psf"
        try:
            import lsst.meas.extensions.psfex.psfexPsfDeterminer
            self.measurePsf.psfDeterminer["psfex"].spatialOrder = 2
            self.measurePsf.psfDeterminer.name = "psfex"
        except ImportError as e:
            print "WARNING: Unable to use psfex: %s" % e
            self.measurePsf.psfDeterminer.name = "pca"

class ReaPsfTask(CmdLineTask):
    ConfigClass = ReaPsfConfig
    _DefaultName = "reaPsf"

    def __init__(self, *args, **kwargs):
        super(ReaPsfTask, self).__init__(*args, **kwargs)
        self.makeSubtask("measurePsf")

    def run(self, dataRef):
        exp = dataRef.get("calexp", immediate=True)
        sources = dataRef.get("icSrc", immediate=True)

        results = self.measurePsf.run(exp, sources)
        psf = results.psf
        cellSet = results.cellSet
        # XXX play with them

        return psf

    def _getConfigName(self):
        pass
    def _getMetadataName(self):
        pass
