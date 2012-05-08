#!/usr/bin/env python

import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

from hsc.pipe.base.argumentParser import HscArgumentParser
from lsst.pex.config import Config, Field, ListField

import numpy
import matplotlib.pyplot as pyplot


class ColorColorConfig(Config):
    filter1 = Field(dtype=str, doc="Name of filter 1", default="W-S-G+")
    filter2 = Field(dtype=str, doc="Name of filter 2", default="W-S-R+")
    filter3 = Field(dtype=str, doc="Name of filter 3", default="W-S-I+")
    filterId = Field(dtype=str, doc="Filter entry in dataId", default="filter")
    fluxes = ListField(dtype=str, doc="List of flux measurements to plot", default=["flux.psf", "flux.sinc"])
    colors = ListField(dtype=str, doc="List of colors", default=["red", "blue", "green", "black"])
    source = Field(dtype=str, doc="Name of data product with sources", default="src")
    radius = Field(dtype=float, doc="Matching radius, arcsec", default=1.0)
    bright1 = Field(dtype=float, doc="Bright limit for filter1", default=22.0)
    bright2 = Field(dtype=float, doc="Bright limit for filter2", default=22.0)
    bright3 = Field(dtype=float, doc="Bright limit for filter3", default=22.0)
    flags = ListField(dtype=str, doc="Flag names to reject",
                      default=["flags.pixel.edge", "flags.pixel.interpolated.center",
                               "flags.pixel.saturated.center"])

    def validate(self):
        if len(self.colors) < len(self.fluxes):
            raise ValueError("Insufficient colors (%d) for fluxes (%d)" % (len(self.colors), len(self.fluxes)))

def newFigure():
    figure = pyplot.figure()

    # Modify pyplot.figure().show() to make it raise the plot too
    def show(self, _show=figure.show):
        _show(self)
        try:
            lift(self)
        except Exception, e:
            pass
            # create a bound method
    import types
    figure.show = types.MethodType(show, figure, figure.__class__)
    return figure

def getCatalog(dataRefList, sourceName, fluxNameList, flags):
    catalogList = [dataRef.get(sourceName) for dataRef in dataRefList]
    zpList = [afwImage.Calib(dataRef.get("calexp_md")).getFluxMag0()[0] for dataRef in dataRefList]
    num = reduce(lambda n, c: n + len(c), catalogList, 0)
    print "Read %d sources" % num

    catalog = afwTable.SourceCatalog(catalogList[0].getSchema())
    catalog.preallocate(num)

    for cat, zp in zip(catalogList, zpList):
        fluxKeyList = [cat.schema.find(fluxName).key for fluxName in fluxNameList]
        flagKeyList = [cat.schema.find(f).key for f in flags]
        classKey = cat.schema.find("classification.extendedness").key
        for record in cat:
            if record.get(classKey) > 0.5: continue
            bad = False
            for key in flagKeyList:
                if record.get(key):
                    bad = True
                    break
            if bad:
                continue

            for key in fluxKeyList:
                record.set(key, record.get(key) / zp)

            catalog.append(catalog.table.copyRecord(record))
    print "Found %d good sources" % len(catalog)      

    return catalog


def colorcolor(dataRefList1, dataRefList2, dataRefList3, config):
    catalog1 = getCatalog(dataRefList1, config.source, config.fluxes, config.flags)
    catalog2 = getCatalog(dataRefList2, config.source, config.fluxes, config.flags)
    catalog3 = getCatalog(dataRefList3, config.source, config.fluxes, config.flags)
    
    match12 = afwTable.matchRaDec(catalog1, catalog2, config.radius * afwGeom.arcseconds)
    match23 = afwTable.matchRaDec(catalog2, catalog3, config.radius * afwGeom.arcseconds)
    print "%d and %d matches" % (len(match12), len(match23))

    sorted12 = sorted(match12, key=lambda m: m.second.getId())
    sorted23 = sorted(match23, key=lambda m: m.first.getId())


    maxLength = max(len(match12), len(match23))

    cat1 = afwTable.SourceCatalog(catalog1.getSchema()); cat1.preallocate(maxLength)
    cat2 = afwTable.SourceCatalog(catalog2.getSchema()); cat2.preallocate(maxLength)
    cat3 = afwTable.SourceCatalog(catalog3.getSchema()); cat3.preallocate(maxLength)

    index23 = 0
    for m1 in sorted12:
        index = m1.second.getId()
        while sorted23[index23].first.getId() < index: index23 += 1
        if sorted23[index23].first.getId() == index:
            cat1.append(cat1.table.copyRecord(m1.first))
            cat2.append(cat2.table.copyRecord(m1.second))
            cat3.append(cat3.table.copyRecord(sorted23[index23].second))
    
    print "Total %d matches" % len(cat1)


    figure = newFigure()
    axes = figure.add_axes((0.1, 0.1, 0.85, 0.80));

    for i, fluxName in enumerate(config.fluxes):
        mag1 = -2.5 * numpy.log10(cat1.columns[fluxName])
        mag2 = -2.5 * numpy.log10(cat2.columns[fluxName])
        mag3 = -2.5 * numpy.log10(cat3.columns[fluxName])

        bright = (mag1 < config.bright1) & (mag2 < config.bright2) & (mag3 < config.bright3)

        axes.plot((mag1-mag2)[bright], (mag2-mag3)[bright], "o", markersize=2,
                  color=config.colors[i], label="%s (%d)" % (fluxName, bright.sum()))
        print "%d good matches for %s" % (bright.sum(), fluxName)

    axes.set_xlabel("%s - %s" % (config.filter1, config.filter2))
    axes.set_ylabel("%s - %s" % (config.filter2, config.filter3))
    axes.legend(loc=2)

    figure.show()



if __name__ == "__main__":
    parser = HscArgumentParser()
    args = parser.parse_args(config=ColorColorConfig())
    
    filterId = args.config.filterId
    dataIdList1 = []
    dataIdList2 = []
    dataIdList3 = []
    for dataRef in args.dataRefList:
        dataId = dataRef.dataId
        filt = dataId[filterId]
        if filt == args.config.filter1:
            dataIdList1.append(dataRef)
            continue
        if filt == args.config.filter2:
            dataIdList2.append(dataRef)
            continue
        if filt == args.config.filter3:
            dataIdList3.append(dataRef)
            continue
        raise RuntimeError("Filter %s is not one of the three recognised filters (%s, %s, %s)" %
                           (filt, args.config.filter1, args.config.filter2, args.config.filter3))

    # Keep plots open when done
    def show():
        print "Please close plots when done."
        try:
            pyplot.show()
        except:
            pass
        print "Plots closed, exiting..."
    import atexit
    atexit.register(show)

    colorcolor(dataIdList1, dataIdList2, dataIdList3, args.config)
