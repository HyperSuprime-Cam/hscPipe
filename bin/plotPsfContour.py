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

def convertImageToNumpyArray(image, exposure, xPsfMapSize=None, yPsfMapSize=None, xGridSize=None, yGridSize=None, xSize=None, ySize=None):
    if xPsfMapSize is None:
        xPsfMapSize = image.getWidth()
    if yPsfMapSize is None:
        yPsfMapSize = image.getHeight()
    if xGridSize is None:
        xGridSize = 1.0
    if yGridSize is None:
        yGridSize = 1.0
    if xSize is None:
        xSize = exposure.getWidth()
    if ySize is None:
        ySize = exposure.getHeight()


    nGridX = int(numpy.floor(float(xSize)/xGridSize))
    nGridY = int(numpy.floor(float(ySize)/yGridSize))
    xPlotSize = xGridSize * nGridX
    yPlotSize = yGridSize * nGridY
    facX = float(xPlotSize)/xPsfMapSize
    facY = float(yPlotSize)/yPsfMapSize

    gridList = numpy.zeros([yPsfMapSize, xPsfMapSize])
    xList, yList = numpy.meshgrid(numpy.linspace(0.5*facX, (xPsfMapSize-0.5)*facX, num=xPsfMapSize, endpoint=True),
                                  numpy.linspace(0.5*facY, (yPsfMapSize-0.5)*facY, num=yPsfMapSize, endpoint=True))
    for j in range(yPsfMapSize):
        for i in range(xPsfMapSize):
            #import pudb; pudb.set_trace()
            #gridList[i, j] = image.get(int(xList[j][i]), int(yList[j][i]))
            gridList[j, i] = image.get(i, j)
            #print gridList

    return  xList, yList, gridList


def getPsfModelGridImage(visit, ccd, psf, exposure, xGridSize=1024, yGridSize=1024, doWriteFits=True, display=True):
    """
    generate an image of psf grid map arranged in a mosaic with given grid sizes
    """
    psfDimension = psf.getKernel().getDimensions()
    xPsfSize = psfDimension.getX()
    yPsfSize = psfDimension.getY()

    #xGridSize = self.config.gridSize
    #yGridSize = self.config.gridSize
    # making grids
    xCcdSize, yCcdSize = exposure.getWidth(), exposure.getHeight()
    nx = int(numpy.floor(float(xCcdSize)/xGridSize))
    ny = int(numpy.floor(float(yCcdSize)/yGridSize))
    print '*** getPsfGridImage: xGridSize, yGridSize:', xGridSize, yGridSize
    print '*** getPsfGridImage: xCcdSize, yCcdSize:', xCcdSize, yCcdSize
    print '*** getPsfGridImage: grid nx, ny:', nx, ny

    # referring to http://dev.lsstcorp.org/doxygen/trunk/afw/classafw_1_1display_1_1utils_1_1_mosaic.html
    m = afwDisp.Mosaic(gutter=0, background=0)
    ###m.setGutter(0) # never call this func; a bug in afwDisp.Mosaic.setGutter() forces to set to 3
    m.setBackground(0)
    #print 'getPsfGridImage: display:nx:', nx

    xc = xGridSize * 0.5
    yc = yGridSize * 0.5
    for j in range(ny):
        for i in range(nx):
            pointXY = afwGeom.Point2D(xc, yc)
            psfSize = afwGeom.Extent2I(xPsfSize, yPsfSize)
            psfImage = psf.computeImage(pointXY, psfSize, True)
            #cellSet.insertCandidate(testSpatialCellLib.ExampleCandidate(xc, yc, psfImage, psfImage.getBBox()))
            m.append(psfImage, '(%d,%d)' % (xc, yc))
            xc += xGridSize
        yc += yGridSize
        xc = xGridSize * 0.5

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

def plotPsfContourGrid(visit, ccd, xList, yList, psfGridList, xPsfMapSize, yPsfMapSize, exposure, xGridSize=1024, yGridSize=1024, xSize=None, ySize=None):
    import matplotlib.figure as figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigCanvas
    if xSize is None:
        xSize = exposure.getWidth()
    if ySize is None:
        ySize = exposure.getHeight()

    facSize = 10.0 / max(xSize,ySize)  # targeting 10 inch in size
    wFig = xSize * facSize * 1.3
    hFig = ySize * facSize
    fig = figure.Figure(figsize=(wFig,hFig))
    canvas = FigCanvas(fig)
    pltPsf = fig.add_axes([0.2, 0.1, 0.7, 0.8]) # left,bottom,width,height

    pltPsf.set_xlim(0, xSize)
    pltPsf.set_ylim(0, ySize)

    #pltPsf.contour(xList, yList, psfGridList)
    #pltPsf.contour(psfGridList, extent=(0, xSize, 0, ySize), origin='lower', extend='neither', levels=[0.01, 0.05,0.25,0.55,0.75,0.95], colors='black', )
    pltPsf.contour(psfGridList, extent=(xList[0][0], xList[-1][-1], yList[0][0], yList[-1][-1]), origin='lower', extend='neither', levels=[0.05,0.25,0.55,0.75,0.95], colors='black', )

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
    parser.add_argument('--xgridsize', type=int, default=1024, help='mesh size in x for psf sampling')
    parser.add_argument('--ygridsize', type=int, default=1024, help='mesh size in y for psf sampling')
    #parser.add_argument('--model', action="store_true", help='use Psf model instead of raw profile')
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


    xGridSize = args.xgridsize
    yGridSize = args.ygridsize

    psfGridImage = getPsfModelGridImage(args.visit, args.ccd, psf, exposure, xGridSize=xGridSize, yGridSize=yGridSize, doWriteFits=True, display=False)

    xPsfMapSize, yPsfMapSize = psfGridImage.getWidth(), psfGridImage.getHeight()
    xCcdSize, yCcdSize = exposure.getWidth(), exposure.getHeight()
    xList, yList, psfGridList = convertImageToNumpyArray(psfGridImage, exposure, xPsfMapSize=xPsfMapSize, yPsfMapSize=yPsfMapSize, xGridSize=xGridSize, yGridSize=yGridSize, xSize=xCcdSize, ySize=yCcdSize)
    #print xList, yList, psfGridList

    plotPsfContourGrid(args.visit, args.ccd, xList, yList, psfGridList, xPsfMapSize, yPsfMapSize, exposure, xGridSize=xGridSize, yGridSize=yGridSize, xSize=xCcdSize, ySize=yCcdSize)

if __name__ == '__main__':
    main()




