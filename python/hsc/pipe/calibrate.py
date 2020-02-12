from lsst.pipe.base import ArgumentParser
from lsst.ctrl.pool.parallel import BatchParallelTask
from lsst.pex.config import Config, Field, ChoiceField
from lsst.afw.table import updateSourceCoords
from lsst.meas.base.forcedPhotCcd import PerTractCcdDataIdContainer
from lsst.afw.geom import calculateSipWcsHeader
from lsst.afw.image import PhotoCalib

__all__ = ("CalibrateCatalogTask", "CalibrateExposureConfig", "CalibrateExposureTask")


class RecalibrateConfig(Config):
    doApplyExternalPhotoCalib = Field(
        dtype=bool,
        default=True,
        doc=("Whether to apply external photometric calibration via an "
             "`lsst.afw.image.PhotoCalib` object. Uses the "
             "``externalPhotoCalibName`` field to determine which calibration "
             "to load."),
    )
    doApplyExternalSkyWcs = Field(
        dtype=bool,
        default=True,
        doc=("Whether to apply external astrometric calibration via an "
             "`lsst.afw.geom.SkyWcs` object. Uses ``externalSkyWcsName`` "
             "field to determine which calibration to load."),
    )
    doApplySkyCorr = Field(
        dtype=bool,
        default=True,
        doc="Apply sky correction?",
    )
    includePhotoCalibVar = Field(
        dtype=bool,
        default=False,
        doc="Add photometric calibration variance to warp variance plane?",
    )
    externalPhotoCalibName = ChoiceField(
        dtype=str,
        doc=("Type of external PhotoCalib if ``doApplyExternalPhotoCalib`` is True. "
             "Unused for Gen3 middleware."),
        default="fgcm_tract",
        allowed={
            "jointcal": "Use jointcal_photoCalib",
            "fgcm": "Use fgcm_photoCalib",
            "fgcm_tract": "Use fgcm_tract_photoCalib"
        },
    )
    externalSkyWcsName = ChoiceField(
        dtype=str,
        doc="Type of external SkyWcs if ``doApplyExternalSkyWcs`` is True. Unused for Gen3 middleware.",
        default="jointcal",
        allowed={
            "jointcal": "Use jointcal_wcs"
        },
    )


class CalibrateCatalogTask(BatchParallelTask):
    ConfigClass = RecalibrateConfig
    _DefaultName = "calibrateCatalog"

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "jointcal_wcs", help="data ID, with raw CCD keys + tract",
                               ContainerClass=PerTractCcdDataIdContainer)
        return parser

    def runDataRef(self, dataRef):
        catalog = dataRef.get("src")

        if self.config.doApplyExternalPhotoCalib:
            source = f"{self.config.externalPhotoCalibName}_photoCalib"
            self.log.info("Applying external photoCalib from %s", source)
            photoCalib = dataRef.get(source)
            catalog = photoCalib.calibrateCatalog(catalog)

        if self.config.doApplyExternalSkyWcs:
            source = f"{self.config.externalSkyWcsName}_wcs"
            self.log.info("Applying external skyWcs from %s", source)
            wcs = dataRef.get(source)
            updateSourceCoords(wcs, catalog)

        dataRef.put(catalog, "calibrated_src")
        return catalog

    def writeConfig(self, *args, **kwargs):
        pass

    def writeSchema(self, *args, **kwargs):
        pass

    def writeMetadata(self, dataRef):
        pass


class CalibrateExposureConfig(RecalibrateConfig):
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

        if self.config.doApplyExternalPhotoCalib:
            source = f"{self.config.externalPhotoCalibName}_photoCalib"
            self.log.info("Applying external photoCalib from %s", source)
            photoCalib = dataRef.get(source)

            exposure.maskedImage = photoCalib.calibrateImage(exposure.maskedImage,
                                                             self.config.includePhotoCalibVar)
            scaling = photoCalib.getCalibrationMean()
            exposure.maskedImage /= scaling
            exposure.setPhotoCalib(PhotoCalib(scaling, photoCalib.getCalibrationErr()))

        if self.config.doApplyExternalSkyWcs:
            source = f"{self.config.externalSkyWcsName}_wcs"
            self.log.info("Applying external skyWcs from %s", source)
            skyWcs = dataRef.get(source)
            exposure.setWcs(skyWcs)
            calculateSipWcsHeader(skyWcs, self.config.sipOrder, exposure.getBBox(), self.config.sipSpacing,
                                  header=exposure.getInfo().getMetadata())

        if self.config.doApplySkyCorr:
            self.log.info("Apply sky correction")
            skyCorr = dataRef.get("skyCorr")
            exposure.maskedImage -= skyCorr.getImage()

        dataRef.put(exposure, "calibrated_exp")
        return exposure

    def writeConfig(self, *args, **kwargs):
        pass

    def writeSchema(self, *args, **kwargs):
        pass

    def writeMetadata(self, dataRef):
        pass

