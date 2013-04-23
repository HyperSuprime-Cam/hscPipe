#!/usr/bin/env python

import numpy
from lsst.pex.config import Config, Field, ListField
from lsst.daf.persistence import Butler
from lsst.afw.image import ImagePcaMF as ImagePca, makeMaskedImage, MaskU, MaskedImageF
import lsst.afw.math as afwMath
from lsst.afw.cameraGeom.utils import ButlerImage, showCamera
import lsst.afw.display.ds9 as ds9
from lsst.pipe.base import ArgumentParser

class PcaConfig(Config):
    """Configuration for PCA"""
    stat = Field(dtype=str, default="MEDIAN", doc="Statistic to use for flux measurement")
    bin = Field(dtype=int, default=4, doc="Binning factor to apply")
    maxMaskedFrac = Field(dtype=float, default=0.1, doc="Maximum tolerable masked fraction")
    bad = ListField(dtype=str, default=["BAD", "SAT",], doc="Mask planes to treat as bad")
    distinguish = Field(dtype=str, default="dateObs", doc="Data identifier name to distinguish nights")

class ButlerMaskedImage(ButlerImage):
    """Interface for showCamera"""
    def __init__(self, *args, **kwargs):
        bg = kwargs.get("background", 0.0)
        super(ButlerMaskedImage, self).__init__(*args, background=bg, **kwargs)
    def getImage(self, ccd, ap=None, imageFactory=None):
        return self.butler.get(self.type, ccd=ccd.getId().getSerial(), **self.kwargs).getMaskedImage()

def addImage(pca, image, bad=["BAD", "SAT"], bgStat="MEDIAN"):
    """Add an image into the PCA"""
    stats = afwMath.StatisticsControl()
    stats.setAndMask(MaskU.getPlaneBitMask(bad))
    flux = afwMath.makeStatistics(image, afwMath.stringToStatisticsProperty(bgStat), stats).getValue()
    pca.addImage(image, flux)

class Null(object):
    """Null object: does nothing"""
    def __init__(self, *args, **kwargs):
        pass
    def __call__(self, *args, **kwargs):
        return self
    def __getattr__(self, name):
        return self
    def __setattr__(self, name, value):
        return self
    def __delattr__(self, name):
        return self

def main(butler, dataRefList, config=PcaConfig(), datasetType="postISRCCD", log=Null()):
    camera = butler.get("camera")
    dataRefList = sorted(dataRefList, cmp=lambda x, y: cmp(x.dataId[config.distinguish],
                                                           y.dataId[config.distinguish]))
    sumImage = None
    lastValue = None
    pca = ImagePca()
    for ref in dataRefList:
        if not ref.datasetExists(datasetType, ccd=0):
            log.warn("Ignoring non-existent data: %s" % ref.dataId)
            continue

        image = showCamera(camera, ButlerMaskedImage(butler, datasetType, **ref.dataId),
                           imageFactory=MaskedImageF, bin=config.bin)
        if (numpy.count_nonzero(numpy.bitwise_and(image.getMask().getArray(),
                                                  MaskU.getPlaneBitMask(config.bad))) > 
            config.maxMaskedFrac * image.getBBox().getArea()):
            log.warn("Ignoring masked image: %s" % ref.dataId)
            continue
        array = image.getImage().getArray()
        indices = numpy.where(numpy.logical_not(numpy.isfinite(array)))
        array[indices] = 0.0

        value = ref.dataId[config.distinguish]
        log.info("Summing: %s" % ref.dataId)
        if sumImage is None:
            sumImage = image
            lastValue = value
        elif value == lastValue:
            sumImage += image
        else:
            #sumImage.writeFits("flat-%03d.fits" % len(pca.getImageList()))
            addImage(pca, sumImage, bgStat=config.stat)
            sumImage = image
            lastValue = value
    addImage(pca, sumImage, bgStat=config.stat)
    del sumImage, image

    num = len(pca.getImageList())
    log.info("Added %d images to PCA" % num)
    pca.analyze()
    pca.updateBadPixels(MaskU.getPlaneBitMask(config.bad), 0)
    pca.analyze()

    eigenValues = pca.getEigenValues()
    eigenImages = pca.getEigenImages()

    sumEigen = sum(eigenValues)
    for i, (image, value) in enumerate(zip(eigenImages, eigenValues)):
        name = "eigen-%03d.fits" % i
        image.writeFits(name)
        log.info("Eigenvalue %d: %e --> %s" % (i, value/sumEigen, name))


if __name__ == "__main__":
    parser = ArgumentParser("flatPca")
    parser.add_id_argument("--id", "raw", level="visit", help="Data identifier")
    parser.add_argument("--datasetType", default="postISRCCD", help="Data set type to read")

    args = parser.parse_args(PcaConfig())

    main(args.butler, args.id.refList, args.config, datasetType=args.datasetType, log=args.log)
