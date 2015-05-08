import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.pipe.base import Task
from lsst.pex.config import Config, ConfigurableField
from lsst.meas.photocal import PhotoCalTask

class PhotometricSolutionConfig(PhotoCalTask.ConfigClass):
    photocal = ConfigurableField(target=PhotoCalTask, doc="Photometric calibration")

    def setDefaults(self):
        self.candidateSourceField = "calib.exposure.photocal.candidate"
        self.usedSourceField = "calib.exposure.photocal.used"

class PhotometricSolutionTask(PhotoCalTask):
    ConfigClass = PhotometricSolutionConfig

    def __init__(self, schema, **kwargs):
        super(PhotometricSolutionTask, self).__init__(schema, **kwargs)
        self.makeSubtask("photocal", schema=schema)

    def run(self, matchLists, filterName):
        matches = self.concatenate(matchLists)
        return self.solve(matches, filterName)

    def concatenate(self, matchLists):
        """Concatenate match lists"""
        matches = []
        for ml in matchLists:
            if ml is None:
                continue
            matches += [m for m in ml]
        return matches

    def solve(self, matches, filterName):
        """Determine a global photometric solution for the exposure.

        The current implementation simply runs the general 'photocal' to get a single zero-point.
        """

        class DummyExposure(object):
            """Quacks like an lsst.afw.image.Exposure, for the purposes of PhotoCal."""
            def __init__(self, filterName):
                self._filterName = filterName
            def getFilter(self):
                return afwImage.Filter(filterName)

        result = self.photocal.run(DummyExposure(filterName), matches)

        self.log.info("Global photometric zero-point: %f" % result.calib.getMagnitude(1.0))
        return result.calib.getFluxMag0()


