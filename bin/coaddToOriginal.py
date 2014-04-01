#!/usr/bin/env python
#!/usr/bin/env python

from argparse import ArgumentParser, Action
from hsc.pipe.tasks.inspectCoadd import coaddToOriginal, identToVisitCcd

class ExposureAction(Action):
    """Supply an ExposureF from a filename"""
    def __call__(self, parser, namespace, filename, option_string):
        import os
        if not os.path.exists(filename):
            parser.error("Exposure %s does not exist" % filename)
        import lsst.afw.image as afwImage
        setattr(namespace, self.dest, afwImage.ExposureF(filename))

class PointAction(Action):
    """Supply a Point2D from x,y"""
    def __call__(self, parser, namespace, values, option_string):
        if len(values) != 2:
            parser.error("Wrong number of coordinates supplied (require 2): %d" % len(values))
        import lsst.afw.geom as afwGeom
        setattr(namespace, self.dest, afwGeom.Point2D(*values))

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("coadd", help="Coadd FITS filename", action=ExposureAction)
    parser.add_argument("point", nargs=2, type=float, help="Coadd coordinates (x,y)", action=PointAction)
    parser.add_argument("--camera", default="hsc", choices=["hsc", "sc"], help="Camera name")
    args = parser.parse_args()

    for ident, point in coaddToOriginal(args.coadd, args.point):
        visit, ccd = identToVisitCcd(args.camera, ident)
        print "visit=%d ccd=%d: %s" % (visit, ccd, point)
