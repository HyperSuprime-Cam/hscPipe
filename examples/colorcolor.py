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
    flux = Field(dtype=str, doc="Name of flux to plot", default="flux.psf")
    source = Field(dtype=str, doc="Name of data product with sources", default="src")
    radius = Field(dtype=float, doc="Matching radius, arcsec", default=1.0)
    bright1 = Field(dtype=float, doc="Bright limit for filter1", default=24.0)
    bright2 = Field(dtype=float, doc="Bright limit for filter2", default=24.0)
    bright3 = Field(dtype=float, doc="Bright limit for filter3", default=24.0)
    flags = ListField(dtype=str, doc="Flag names to reject",
                      default=["flags.pixel.edge", "flags.pixel.interpolated.center",
                               "flags.pixel.saturated.center"])

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

def getCatalog(dataRefList, sourceName, fluxName, flags):
    catalogList = [dataRef.get(sourceName) for dataRef in dataRefList]
    zpList = [afwImage.Calib(dataRef.get("calexp_md")).getFluxMag0()[0] for dataRef in dataRefList]
    num = reduce(lambda n, c: n + len(c), catalogList, 0)
    print "Read %d sources" % num

    catalog = afwTable.SourceCatalog(catalogList[0].getSchema())
    catalog.preallocate(num)

    for cat, zp in zip(catalogList, zpList):
        fluxKey = cat.schema.find(fluxName).key
        flagKeys = [cat.schema.find(f).key for f in flags]
        classKey = cat.schema.find("classification.extendedness").key
        for record in cat:
            if record.get(classKey) > 0.5: continue
            bad = False
            for k in flagKeys:
                if record.get(k):
                    bad = True
                    break
            if bad:
                continue

            record.set(fluxKey, record.get(fluxKey) / zp)
            catalog.append(catalog.table.copyRecord(record))
    print "Found %d good sources" % len(catalog)      

    return catalog


def colorcolor(dataRefList1, dataRefList2, dataRefList3, config):
    catalog1 = getCatalog(dataRefList1, config.source, config.flux, config.flags)
    catalog2 = getCatalog(dataRefList2, config.source, config.flux, config.flags)
    catalog3 = getCatalog(dataRefList3, config.source, config.flux, config.flags)
    
    match12 = afwTable.matchRaDec(catalog1, catalog2, config.radius * afwGeom.arcseconds)
    match23 = afwTable.matchRaDec(catalog2, catalog3, config.radius * afwGeom.arcseconds)
    print "%d and %d matches" % (len(match12), len(match23))

    sorted12 = sorted(match12, key=lambda m: m.second.getId())
    sorted23 = sorted(match12, key=lambda m: m.first.getId())


    maxLength = max(len(match12), len(match23))

    cat1 = afwTable.SourceCatalog(catalog1.getSchema()); cat1.preallocate(maxLength)
    cat2 = afwTable.SourceCatalog(catalog2.getSchema()); cat2.preallocate(maxLength)
    cat3 = afwTable.SourceCatalog(catalog3.getSchema()); cat3.preallocate(maxLength)

    index23 = 0
    for m1 in sorted12:
        index = m1.second.getId()
        while match23[index23].first.getId() < index: index23 += 1
        if match23[index23].first.getId() == index:
            cat1.append(cat1.table.copyRecord(m1.first))
            cat2.append(cat2.table.copyRecord(m1.second))
            cat3.append(cat3.table.copyRecord(match23[index23].second))
    
    print "Total %d matches" % len(cat1)

    mag1 = -2.5 * numpy.log10(cat1.columns[config.flux])
    mag2 = -2.5 * numpy.log10(cat2.columns[config.flux])
    mag3 = -2.5 * numpy.log10(cat3.columns[config.flux])

    bright = (mag1 < config.bright1) & (mag2 < config.bright2) & (mag3 < config.bright3)

    figure = newFigure()
    axes = figure.add_axes((0.1, 0.1, 0.85, 0.80));
    axes.plot((mag1-mag2)[bright], (mag2-mag3)[bright], "o", markersize=2, color="red")

    axes.set_xlabel("%s - %s" % (config.filter1, config.filter2))
    axes.set_ylabel("%s - %s" % (config.filter2, config.filter3))
    axes.set_title(config.flux)

    fig.show()



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

    colorcolor(dataIdList1, dataIdList2, dataIdList3, args.config)
