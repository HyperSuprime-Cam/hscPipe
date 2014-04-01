import functools
from lsst.pex.config import Config, Field
from lsst.pipe.base import CmdLineTask, ArgumentParser
from lsst.pipe.tasks.coaddBase import ExistingCoaddDataIdContainer
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9

def coaddToOriginal(exposure, point):
    """Convert coadd coordinates to original coordinates"""
    info = exposure.getInfo()
    if not info:
        raise RuntimeError("No exposure info")
    inputs = info.getCoaddInputs()
    if not inputs:
        raise RuntimeError("No coadd inputs")
    wcs = exposure.getWcs()
    if not wcs:
        raise RuntimeError("No coadd Wcs")
    sky = wcs.pixelToSky(point)
    return [(orig.getId(), orig.getWcs().skyToPixel(sky)) for orig in inputs.ccds if orig.contains(sky)]

def identToVisitCcd(camera, ident):
    """Convert identifier to visit,ccd pair"""
    if camera.lower() == "hsc":
        return ident//200, ident % 200
    if camera.lower() == "sc":
        return ident//10, ident % 10
    raise RuntimeError("Unrecognised camera name: %s" % camera)


class InspectCoaddConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd type (without the 'Coadd' suffix)")
    origName = Field(dtype=str, default="calexp", doc="Name of original type")
    maskTransparency = Field(dtype=float, default=90, doc="ds9 mask transparency")
    zoom = Field(dtype=float, default=4.0, doc="Zoom to use on originals")

class InspectCoaddTask(CmdLineTask):
    ConfigClass = InspectCoaddConfig
    _DefaultName = "inspect"

    @classmethod
    def _makeArgumentParser(cls):
        parser = ArgumentParser(name=cls._DefaultName)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=ExistingCoaddDataIdContainer)
        return parser

    def run(self, dataRef):
        self.butler = dataRef.getButler()
        self.frames = 0

        exposure = dataRef.get()
        ds9.setCallback(" ", functools.partial(self.inspect, exposure))
        ds9.setCallback("i", functools.partial(self.inspect, exposure))
        ds9.setCallback("q", lambda k, x, y: True)

        # Suppress annoying "No callback is registered for <key>" messages
        nullCallback = lambda k, x, y: False
        ds9.setCallback("Tab", nullCallback)
        ds9.setCallback("Shift_L", nullCallback)
        ds9.setCallback("ISO_Left_tab", nullCallback)
        ds9.setCallback("ISO_Left_Tab", nullCallback)

        ds9.setMaskTransparency(self.config.maskTransparency)
        ds9.mtv(exposure, title="Coadd", frame=1)
        self.log.info("ds9 display is interactive: 'i' to inspect position on original images, 'q' to quit")
        ds9.interact()

    def inspect(self, exposure, key, x, y):
        originalList = coaddToOriginal(exposure, afwGeom.Point2D(x, y) + afwGeom.Extent2D(exposure.getXY0()))
        if len(originalList) == 0:
            self.log.info("No originals at coadd point %f,%f" % (x,y))
            return False
        self.log.info("Displaying originals at coadd point %f,%f" % (x,y))
        for i, (ident, point) in enumerate(originalList):
            visit, ccd = identToVisitCcd(self.butler.mapper.getCameraName(), ident)
            image = self.butler.get(self.config.origName, visit=visit, ccd=ccd)
            frame = i + 2
            self.log.info("Frame %d: visit=%d ccd=%d: %s" % (frame, visit, ccd, point))
            ds9.mtv(image, title="visit=%d ccd=%d" % (visit, ccd), frame=frame)
            point -= afwGeom.Extent2D(image.getXY0())
            ds9.dot("o", point.getX(), point.getY(), frame=frame)
            ds9.zoom(self.config.zoom, point.getX(), point.getY(), frame=frame)
        if self.frames > len(originalList):
            for i in range(len(originalList), self.frames):
                ds9.ds9Cmd("frame delete %d" % (i + 2))
        self.frames = len(originalList)
        self.log.info("Return to coadd (frame=1) to continue inspecting")
        return False

    def _getConfigName(self):
        return None
    def _getMetadataName(self):
        return None

