#!/usr/bin/env python

import os, sys
import hsc.pipe.base.camera  as hscCamera
import lsst.meas.algorithms as measAlg
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.display.ds9 as ds9
import lsst.afw.display.utils as afwDisp

import numpy


def convertImageToNumpyArray2(image, x0=0, y0=0, xmax=None, ymax=None):
    if xmax is None:
        xmax = image.getWidth()
    if ymax is None:
        ymax = image.getHeight()
    gridArray = numpy.zeros((xmax, ymax))
    imbb = image.getBBox(afwImage.PARENT)
    imx0, imy0 = imbb.getMinX(), imbb.getMinY()
    for j in range(ymax):
        for i in range(xmax):
            if imbb.contains(afwGeom.Point2I(i+x0, j+y0)):
                gridArray[i, j] = image.get(i+x0-imx0, j+y0-imy0)
    return gridArray

def convertImageToNumpyArray(image, nx=None, ny=None):
    if nx is None:
        nx = image.getWidth()
    if ny is None:
        ny = image.getHeight()
    gridList = numpy.zeros([ny, nx])
    xList, yList = numpy.meshgrid(numpy.arange(0, nx, 1), numpy.arange(0, ny, 1))
    for j in range(ny):
        for i in range(nx):
            #import pudb; pudb.set_trace()
            #gridList[i, j] = image.get(int(xList[j][i]), int(yList[j][i]))
            gridList[j, i] = image.get(i, j)
            print gridList

    return  xList, yList, gridList

if False:
    def plotPsfContourGrid(self, dataRef, psf, exposure):
        """fwhm-in-grids plots"""
        xSize, ySize = exposure.getWidth(), exposure.getHeight()
        facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
        wFig = xSize * facSize * 1.3
        hFig = ySize * facSize
        fig = figure.Figure(figsize=(wFig,hFig))
        canvas = FigCanvas(fig)
        pltFwhm = fig.add_axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height

        pltFwhm.set_xlim(0, xSize)
        pltFwhm.set_ylim(0, ySize)

        xGridSize = self.config.gridSize
        yGridSize = self.config.gridSize

        # making grids
        nx = int(numpy.floor(float(xSize)/xGridSize))
        ny = int(numpy.floor(float(ySize)/yGridSize))

        xGrids = numpy.array([ (i+0.5)*xGridSize for i in range(nx) ])
        yGrids = numpy.array([ (i+0.5)*yGridSize for i in range(ny) ])
        xMeshGrid, yMeshGrid = numpy.meshgrid(xGrids, yGrids)
        xGridList = numpy.reshape(xMeshGrid, (-1,))
        yGridList = numpy.reshape(yMeshGrid, (-1,))

        # walking through all grid meshes
        fwhmList = numpy.array([])
        for xGridCenter, yGridCenter in zip(xGridList, yGridList):
            xGridMin = xGridCenter - xGridSize/2.
            xGridMax = xGridCenter + xGridSize/2.
            yGridMin = yGridCenter - yGridSize/2.
            yGridMax = yGridCenter + yGridSize/2.

            fwhmsInGrid = numpy.array([])
            for xFwhm, yFwhm, fwhm in zip(data.xListPsfLikeRobust, data.yListPsfLikeRobust,
                                          data.fwhmListPsfLikeRobust):
                if xGridMin <= xFwhm and xFwhm < xGridMax and yGridMin <= yFwhm and yFwhm < yGridMax:
                    if fwhm is not None:
                        fwhmsInGrid = numpy.append(fwhmsInGrid, fwhm)
            # taking median to represent the value in a grid mesh
            fwhmList = numpy.append(fwhmList, numpy.median(fwhmsInGrid))

        # For latter use; this may be a bit entangled way.
        data.nxGrid = nx
        data.nyGrid = ny
        data.yGridList = yGridList
        data.xGridList = xGridList
        data.yGridList = yGridList
        data.fwhmGridList = fwhmList

        # 10pix=2arcsec(fwhm)=500 point(radius) (to be ~0.6*min(xGridSize, yGridSize)?)
        pointRadius = 100*fwhmList/2.
        scaleFactor = min(xGridSize/xSize, yGridSize/ySize)
        pointRadius *= scaleFactor    
        pointArea = math.pi*(pointRadius)**2.

        pltFwhm.scatter(xGridList, yGridList, s=pointArea, marker='o', color=None, facecolor=(1,1,1,0),
                        linewidth=5.0, label='PSF sample')

        # reference sample symbol
        fwhmPix = 5.0 # 5pix in fwhm
        pointRadius = 100*numpy.array([fwhmPix])/2. 
        scaleFactor = min(xGridSize/xSize, yGridSize/ySize)
        pointRadius *= scaleFactor    
        pointArea = math.pi*(pointRadius)**2.
        pltFwhm.scatter([0.2 * xSize], [0.85 * ySize], s=pointArea, marker='o', color='magenta',
                        facecolor=(1,1,1,0), linewidth=8.0, label='PSF sample')
        fig.text(0.2 * xSize, 0.9 * ySize, 'fwhm=%4.1f pix' % fwhmPix, ha='center', va='top')

        pltFwhm.set_title('FWHM of PSF sources')
        pltFwhm.set_xlabel('X (pix)')
        pltFwhm.set_ylabel('Y (pix)')

        pltFwhm.set_xticks([ xc+xGridSize/2. for xc in xGridList ])
        pltFwhm.set_yticks([ yc+yGridSize/2. for yc in yGridList ])
        pltFwhm.grid()

        fname = getFilename(dataRef, "plotFwhmGrid")
        fig.savefig(fname, dpi=None, facecolor='w', edgecolor='w', orientation='portrait',
                    papertype=None, format='png', transparent=False, bbox_inches=None, pad_inches=0.1)

        del fig
        del pltFwhm
        del canvas

def getPsfGridImage(visit, ccd, psf, exposure, xGridSize=1024, yGridSize=1024, doWriteFits=True, display=True):
    """
    generate an image of psf grid map arranged in a mosaic with given grid sizes
    """
    psfDimension = psf.getKernel().getDimensions()
    xSizePsf = psfDimension.getX()
    ySizePsf = psfDimension.getY()

    #xGridSize = self.config.gridSize
    #yGridSize = self.config.gridSize
    # making grids
    xSize, ySize = exposure.getWidth(), exposure.getHeight()
    print '** xSize, ySize:', xSize, ySize
    nx = int(numpy.floor(float(xSize)/xGridSize))
    ny = int(numpy.floor(float(ySize)/yGridSize))

    # referring to http://dev.lsstcorp.org/doxygen/trunk/afw/classafw_1_1display_1_1utils_1_1_mosaic.html
    if True: # if config.display:
        m = afwDisp.Mosaic(gutter=0, background=0)
        ###m.setGutter(0) # never call this func; a bug in afwDisp.Mosaic.setGutter() forces to set to 3
        m.setBackground(0)
        print 'display:nx:', nx

    xc = xGridSize * 0.5
    yc = yGridSize * 0.5
    while yc < ySize:
        while xc < xSize:
            pointXY = afwGeom.Point2D(xc, yc)
            psfSize = afwGeom.Extent2I(xSizePsf, ySizePsf)
            psfImage = psf.computeImage(pointXY, psfSize, True)

            xc += xGridSize
            #cellSet.insertCandidate(testSpatialCellLib.ExampleCandidate(xc, yc, psfImage, psfImage.getBBox()))
            if True: #if config.display:
                m.append(psfImage, '(%d,%d)' % (xc, yc))

        xc = xGridSize * 0.5
        yc += yGridSize

    m.drawLabels()
        # See afw.display.utils.py:L129.
        #    Mosaic.makeMoasic() seems to accept nx*ny mosaic, in addition to 'square', 'x', and 'y'.
        #    m.setMode(nx)
    psfGridImage = m.makeMosaic(mode=nx)

    if display:
        ds9.mtv(psfGridImage, frame=0)

    if doWriteFits:
        psfGridImage.writeFits('psfGrid-%07d-%03d.fits' % (visit, ccd))

    return psfGridImage


def plotPsfCountorGrid(visit, ccd, xList, yList, psfGridList, xSize, ySize, exposure, xGridSize=1024, yGridSize=1024):
    import matplotlib.figure as figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas

    #xSize, ySize = exposure.getWidth(), exposure.getHeight()
    facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
    wFig = xSize * facSize * 1.3
    hFig = ySize * facSize
    fig = figure.Figure(figsize=(wFig,hFig))
    canvas = FigCanvas(fig)
    pltPsf = fig.add_axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height

    pltPsf.set_xlim(0, xSize)
    pltPsf.set_ylim(0, ySize)

    #pltPsf.contour(xList, yList, psfGridList)
    pltPsf.contour(psfGridList, extent=(0, xSize, 0, ySize), origin='lower', extend='neither', levels=[0.01, 0.05,0.25,0.55,0.75,0.95], colors='black', )

    pltPsf.set_title('Psf Contour')
    pltPsf.set_xlabel('X (pix)')
    pltPsf.set_ylabel('Y (pix)')

    filename = 'psfGrid-%07d-%03d.png' % (visit, ccd)
    fig.savefig(filename, dpi=None, facecolor='w', edgecolor='w', orientation='portrait', papertype=None,
                format='png', transparent=False, bbox_inches=None, pad_inches=0.1)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('instrument')
    parser.add_argument('rerun')
    parser.add_argument('visit', type=int)
    parser.add_argument('ccd', type=int)
    parser.add_argument('--inroot', default='/data/data1', help='e.g., /data/data1')
    parser.add_argument('--outroot', default='/data/data2', help='e.g., /data/data2/')
    #parser.add_argument('--dst-rerun', dest='dstRerun', default=None, help='XXX')
    args = parser.parse_args()

    if args.instrument.lower() in ['hsc', 'hscsim']:
        addDir = 'HSC'
        if args.instrument.lower() == 'hsc' or args.instrument == 'hscsim':
            args.instrument = 'hscSim'
    else:
        addDir = 'SUPA'

    inRootDir = os.path.join(args.inroot, 'Subaru', addDir)
    outRootDir = os.path.join(args.outroot, 'Subaru', addDir, 'rerun', args.rerun)
    calibRootDir = os.path.join(inRootDir, 'CALIB')

    butler = hscCamera.getButler(args.instrument, rerun=args.rerun, root=inRootDir, outputRoot=outRootDir, calibRoot=calibRootDir)
    exposure = butler.get('calexp', {'visit':args.visit, 'ccd':args.ccd})
    psf = butler.get('psf', {'visit':args.visit, 'ccd':args.ccd})


    gridSize = 2048

    psfGridImage = getPsfGridImage(args.visit, args.ccd, psf, exposure, xGridSize=gridSize, yGridSize=gridSize, doWriteFits=True, display=False)


    xPsfMapSize, yPsfMapSize = psfGridImage.getWidth(), psfGridImage.getHeight()
    xList, yList, psfGridList = convertImageToNumpyArray(psfGridImage, nx=xPsfMapSize, ny=yPsfMapSize)
    print xList, yList, psfGridList

    plotPsfCountorGrid(args.visit, args.ccd, xList, yList, psfGridList, xPsfMapSize, yPsfMapSize, exposure, xGridSize=gridSize, yGridSize=gridSize)


if __name__ == '__main__':
    main()




