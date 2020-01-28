from lsst.pipe.base import ArgumentParser
from lsst.ctrl.pool.parallel import BatchParallelTask
from lsst.pex.config import Config, Field
from lsst.afw.table import updateSourceCoords
from lsst.meas.base.forcedPhotCcd import PerTractCcdDataIdContainer
from lsst.afw.geom import calculateSipWcsHeader

__all__ = ("CalibrateCatalogTask", "CalibrateExposureConfig", "CalibrateExposureTask")


class CalibrateCatalogTask(BatchParallelTask):
    ConfigClass = Config
    _DefaultName = "calibrateCatalog"

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "jointcal_wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    def runDataRef(self, dataRef):
        catalog = dataRef.get("src")

        photoCalib = dataRef.get("jointcal_photoCalib")
        catalog = photoCalib.calibrateCatalog(catalog)

        wcs = dataRef.get("jointcal_wcs")
        updateSourceCoords(wcs, catalog)

        dataRef.put(catalog, "calibrated_src")
        return catalog

    def writeConfig(self, *args, **kwargs):
        pass

    def writeSchema(self, *args, **kwargs):
        pass

    def writeMetadata(self, dataRef):
        pass


class CalibrateExposureConfig(Config):
    includeCalibVar = Field(dtype=bool, default=False, doc="Add photometric calibration variance?")
    sipOrder = Field(dtype=int, default=9, doc="Polynomial order for SIP WCS")
    sipSpacing = Field(dtype=float, default=20, doc="Spacing in pixels between samples for SIP fit")


class CalibrateExposureTask(BatchParallelTask):
    ConfigClass = CalibrateExposureConfig
    _DefaultName = "calibrateExposure"

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "jointcal_wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    def runDataRef(self, dataRef):
        exposure = dataRef.get("calexp")

        photoCalib = dataRef.get("jointcal_photoCalib")
        exposure.maskedImage = photoCalib.calibrateImage(exposure.maskedImage, self.config.includeCalibVar)

        wcs = dataRef.get("jointcal_wcs")
        exposure.setWcs(wcs)
        calculateSipWcsHeader(wcs, self.config.sipOrder, exposure.getBBox(), self.config.sipSpacing,
                              header=exposure.getInfo().getMetadata())

        dataRef.put(exposure, "calibrated_exp")
        return exposure

    def writeConfig(self, *args, **kwargs):
        pass

    def writeSchema(self, *args, **kwargs):
        pass

    def writeMetadata(self, dataRef):
        pass

