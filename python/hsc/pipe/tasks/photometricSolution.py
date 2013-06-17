import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.pipe.base import Task
from lsst.pex.config import Config, ConfigurableField
from lsst.meas.photocal import PhotoCalTask

class PhotometricSolutionConfig(PhotoCalTask.ConfigClass):
    photocal = ConfigurableField(target=PhotoCalTask, doc="Photometric calibration")

    def setDefaults(self):
        self.outputField = "classification.exposure.photometric"

class PhotometricSolutionTask(PhotoCalTask):
    ConfigClass = PhotometricSolutionConfig

    def __init__(self, schema, **kwargs):
        super(PhotometricSolutionTask, self).__init__(schema, **kwargs)
        self.makeSubtask("photocal", schema=schema)

    def run(self, matchLists, filterName):
        matches = self.concatenate(matchLists)
        return self.solve(matches, filterName)

    def concatenate(self, matchLists):
        """Concatenate match lists

        We take some care that the schemas differ.  It's possible
        the SchemaMapper would help here, but this works for now.
        """
        template = None
        for ml in matchLists:
            if ml is not None and len(ml) > 0:
                template = ml[0]
                break
        if template is None:
            self.log.warn("No matches provided; setting crazy zero point")
            return 0.0
        ref = afwTable.SimpleTable.make(template.first.schema)
        src = afwTable.SourceTable.make(template.second.schema)
        for prop in ("Centroid", "Shape", "PsfFlux", "ApFlux", "ModelFlux", "InstFlux"):
            getter = getattr(template.second.table, "get" + prop + "Definition")
            setter = getattr(src, "define" + prop)
            setter(getter())

        refNames = ref.schema.getNames()
        srcNames = src.schema.getNames()
        refNewKeys = [ref.schema.find(k).key for k in refNames]
        srcNewKeys = [src.schema.find(k).key for k in srcNames]
        matches = []
        for i, ml in enumerate(matchLists):
            if ml is None or len(ml) == 0:
                continue
            try:
                refOldKeys = [ml[0].first.schema.find(k).key for k in refNames]
                srcOldKeys = [ml[0].second.schema.find(k).key for k in srcNames]
            except:
                # Something's wrong with the schema; punt
                self.log.warn("Error with schema on matchlist %d: ignoring %d matches\n" % (i, len(ml)))
                continue

            for m in ml:
                newRef = ref.makeRecord()
                for old, new in zip(refOldKeys, refNewKeys):
                    newRef.set(new, m.first.get(old))
                newSrc = src.makeRecord()
                for old, new in zip(srcOldKeys, srcNewKeys):
                    newSrc.set(new, m.second.get(old))
                matches.append(afwTable.ReferenceMatch(newRef, newSrc, m.distance))
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


