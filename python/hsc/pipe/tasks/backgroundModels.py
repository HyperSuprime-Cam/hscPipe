import numpy
from lsst.pex.config import Config, Field, ListField, ChoiceField
import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath


BINNING = 20

def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.mean(%d)'%(i+1) for i in range(lenShape)]
    return eval(''.join(evList))

class Fitter(object):
    # Least squares fit:
    # * Matrix contains sum(model.dot.model)
    # * Vector contains sum(model.dot.data)
    # Note that even though we're using chebyshevs, model.dot.model is not zero because we're not
    # integrating over the entire range (e.g., we have a mask).

    def __init__(self, numParams):
        self.numParams = numParams
        self.matrix = numpy.zeros((numParams, numParams), dtype=numpy.float64)
        self.vector = numpy.zeros((numParams), dtype=numpy.float64)

    def accumulate(self, image, evaluator, badPixels=["BAD", "SAT", "INTRP", "DETECTED", "EDGE",]):
        box = image.getBBox(afwImage.PARENT)
        xMin, yMin = box.getMin()
        xMax, yMax = box.getMax()

        if hasattr(image, "getMask"):
            maskVal = image.getMask().getPlaneBitMask(badPixels)
            mask = (image.getMask().getArray() & maskVal > 0)
            data = numpy.ma.masked_where(mask, image.getImage().getArray())
        else:
            data = image.getArray()
        data = numpy.ma.masked_where(numpy.logical_not(numpy.isfinite(data)), data)
        print data.mean()
        data = rebin(data, data.shape[0]/BINNING, data.shape[1]/BINNING)
        print data.mean()

        for iParam in range(self.numParams):
            iFunc = evaluator[iParam]
            self.vector[iParam] += (iFunc*data).sum()
            self.matrix[iParam,iParam] += (iFunc*iFunc).sum()
            for jParam in range(iParam + 1, self.numParams):
                jFunc = evaluator[jParam]
                ijValue = (iFunc*jFunc).sum()
                self.matrix[iParam,jParam] += ijValue
                self.matrix[jParam,iParam] = self.matrix[iParam,jParam]
            print "%d/%d" % (iParam, self.numParams)
        return self

    def finish(self):
        solution = afwMath.LeastSquares.fromNormalEquations(self.matrix, self.vector).getSolution()
        return solution.astype(numpy.float32)


class FunctionEvaluator(object):
    def __init__(self, xFunc, yFunc, box, binning=BINNING):
        self.box = box
        xNum, yNum = xFunc.getNParameters(), yFunc.getNParameters()
        self.numParams = xNum*yNum
        self._xIndices = numpy.ndarray(self.numParams, dtype=int)
        self._yIndices = numpy.ndarray(self.numParams, dtype=int)
        i = 0
        for y in range(yNum):
            for x in range(xNum):
                self._xIndices[i] = x
                self._yIndices[i] = y
                i += 1

        def evaluateSingle(func, low, high):
            num = func.getNParameters()
            size = high - low + 1
            params = numpy.zeros(num)
            evaluation = numpy.ndarray((num, size/binning))
            values = numpy.arange(low + 0.5*(binning - 1), high + 1, binning)
            for i in range(num):
                params[i] = 1.0
                func.setParameters(params)
                evaluation[i,:] = [func(x) for x in values] # operator() isn't vectorised!
                params[i] = 0.0
            return evaluation

        self._x = evaluateSingle(xFunc, box.getMinX(), box.getMaxX())
        self._y = evaluateSingle(yFunc, box.getMinY(), box.getMaxY())

    def __getitem__(self, iParam):
        """Get evaluation for just this parameter"""
        return numpy.outer(self._y[self._yIndices[iParam]], self._x[self._xIndices[iParam]])

    def __call__(self, coeffs):
        """Evaluate with these coefficients"""
        return reduce(lambda image, i: image + coeffs[i]*self[i], range(1, self.numParams), coeffs[0]*self[0])


class ChebyBackground(object):
    def __init__(self, xOrder, yOrder, box, coeffs=None):
        Function = afwMath.Chebyshev1Function1F
        self.xOrder = xOrder
        self.yOrder = yOrder
        self.box = box
        xMin, yMin = box.getMin()
        xMax, yMax = box.getMax()
        self._xFunc = Function(xOrder, xMin, xMax)
        self._yFunc = Function(yOrder, yMin, yMax)
        self.numParams = self._xFunc.getNParameters()*self._yFunc.getNParameters()
        if coeffs is not None:
            if len(coeffs) != self.numParams:
                raise RuntimeError("Number of coefficients doens't match expected: %d vs %d" %
                                   (len(coeffs), self.numParams))
        else:
            coeffs = numpy.zeros(self.numParams)
        self.coeffs = coeffs

    def __reduce__(self):
        return self.__class__, (self.xOrder, self.yOrder, self.box, self.coeffs)

    @classmethod
    def fromImage(cls, image, xOrder, yOrder, badPixels=["BAD", "SAT", "INTRP", "DETECTED", "EDGE",], box=None):
        if hasattr(image, "getMaskedImage"):
            image = image.getMaskedImage()
        if box is None:
            box = image.getBBox(afwImage.PARENT)
        self = cls(xOrder, yOrder, box)
        evaluator = self._getEvaluator(image.getBBox(afwImage.PARENT))
        self.coeffs = Fitter(self.numParams).accumulate(image, evaluator, badPixels).finish()
        print "Coefficients: %s" % self.coeffs
        return self

    def _getEvaluator(self, box=None, binning=BINNING):
        if box is None:
            box = self.box
        return FunctionEvaluator(self._xFunc, self._yFunc, box, binning=binning)

    def getImage(self, *args):
        # Reproduce API for BackgroundMI.getImageF
        box = None
        for a in args:
            if isinstance(a, afwGeom.Box2I):
                box = a

        if box is None:
            box = self.box
        image = afwImage.ImageF(box)
        image.getArray()[:] = self._getEvaluator(box, binning=1)(self.coeffs)
        return image

    getImageF = getImage

    @classmethod
    def fromMerge(cls, bgList, xOrder, yOrder, box=None):
        if box is None:
            box = afwGeom.Box2I()
            for bg in bgList:
                box.include(bg.box)

        self = cls(xOrder, yOrder, box)
        fitter = Fitter(self.numParams)
        for i, bg in enumerate(bgList):
            print "Accumulating %d: %s" % (i, bg.box)
            image = bg.getImage()
            fitter.accumulate(image, self._getEvaluator(bg.box))
        self.coeffs = fitter.finish()
        print "Coefficients: %s" % self.coeffs
        return self


def getBoundsCenters(low, high, size):
    # Bounds are inclusive (i.e., suitable for Box2I)
    # Over-sampling
    numBounds = 2*(high - low + 1)//size + 1
    posts = numpy.rint(numpy.linspace(low, high + 1, num=numBounds, endpoint=True)).astype(int)
    bounds = zip(posts[:-2], posts[2:] - 1)
    return bounds, numpy.array([0.5*(b2 + b1) for b1, b2 in bounds])

class BackgroundConfig(Config):
    xSize = Field(dtype=int, default=256, doc="Binning size in x")
    ySize = Field(dtype=int, default=256, doc="Binning size in y")
    minFrac = Field(dtype=float, default=0.1, doc="Minimum fraction of bin size for good measurement")
    badPixels = ListField(dtype=str, default=["BAD", "SAT", "INTRP", "DETECTED", "EDGE",],
                          doc="Mask planes to treat as bad")
    maxExtrapolateFrac = Field(dtype=float, default=0.5,
                               doc="Maximum fraction of a bin to extrapolate; 0 for unbounded")
    extrapolateValue = Field(dtype=float, default=numpy.nan, doc="Value for extrapolation if disallowed")
    extrapolationMode = ChoiceField(dtype=str, doc="Mode of extrapolation", default="NONE",
                                    allowed={"NONE": "No extrapolation allowed",
                                             "CONSTANT": "Extrapolate with a constant",
                                             "LIMITED": "Limit extrapolation with maxExtrapolateFrac",
                                             "FULL": "Full extrapolation allowed",
                                             }, optional=False)
    approximateStyle = ChoiceField(dtype=str, doc="Approximation style for statistics image", default="NONE",
                                   allowed={"NONE": "No approximation",
                                            "CHEBYSHEV": "Approximate with Chebyshev polynomial",
                                        }, optional=False)
    xApproximate = Field(dtype=int, default=5, doc="Approximation order in x")
    yApproximate = Field(dtype=int, default=5, doc="Approximation order in y")

class Background(object):
    def __init__(self, config, box, polygon=None, values=None, numbers=None):
        self.box = box
        self.config = config
        self.polygon = polygon

        self._xBounds, self._xCenter = getBoundsCenters(box.getMinX(), box.getMaxX(), config.xSize)
        self._yBounds, self._yCenter = getBoundsCenters(box.getMinY(), box.getMaxY(), config.ySize)

        xNum, yNum = len(self._xCenter), len(self._yCenter)

        if False:
            print self._xBounds
            print self._xCenter
            print self._yBounds
            print self._yCenter
            for (yMin, yMax), yCenter in zip(self._yBounds, self._yCenter):
                for (xMin, xMax), xCenter in zip(self._xBounds, self._xCenter):
                    b = afwGeom.Box2I(afwGeom.Point2I(xMin, yMin), afwGeom.Point2I(xMax, yMax))
                    print b, b.getDimensions(), afwGeom.Box2D(b).getCenter(), afwGeom.Point2D(xCenter, yCenter)

        imageBox = afwGeom.Box2I(afwGeom.Point2I(0,0),
                                 afwGeom.Point2I(len(self._xCenter) - 1, len(self._yCenter) - 1)) # inclusive
        if values is None:
            values = afwImage.ImageF(imageBox)
            values.set(0.0)
        assert(values.getBBox(afwImage.PARENT) == imageBox)
        assert(values.getWidth() == xNum and values.getHeight() == yNum)
        self._values = values
        if numbers is None:
            numbers = afwImage.ImageF(imageBox) # float for dynamic range and convenience
            numbers.set(0.0)
        assert(numbers.getBBox(afwImage.PARENT) == imageBox)
        self._numbers = numbers

    def __reduce__(self):
        return self.__class__, (self.config, self.box, self.polygon, self._values, self._numbers)

    def addImage(self, image):
        if hasattr(image, "getMaskedImage"):
            image = image.getMaskedImage()
        imageBox = image.getBBox(afwImage.PARENT)
        stats = afwMath.StatisticsControl()
        statistic = afwMath.MEAN
        mask = None
        if hasattr(image, "getMask"):
            mask = image.getMask()
            maskVal = mask.getPlaneBitMask(self.config.badPixels)
            stats.setAndMask(maskVal)
        stats.setNanSafe(True)

        if self.polygon is not None:
            polyImage = self.polygon.createImage(imageBox)

        xMinImage, yMinImage = imageBox.getMin()
        xMaxImage, yMaxImage = imageBox.getMax()

        for j, (yMin, yMax) in enumerate(self._yBounds):
            if yMax < yMinImage or yMin > yMaxImage:
                continue
            for i, (xMin, xMax) in enumerate(self._xBounds):
                if xMax < xMinImage or xMin > xMaxImage:
                    continue
                box = afwGeom.Box2I(afwGeom.Point2I(int(xMin), int(yMin)),
                                    afwGeom.Point2I(int(xMax), int(yMax)))
                box.clip(imageBox)
                if self.polygon is not None and not self.polygon.overlaps(afwGeom.Box2D(box)):
                    continue
                subImage = image.Factory(image, box, afwImage.PARENT)
                if self.polygon is not None:
                    subPoly = polyImage.Factory(polyImage, box, afwImage.PARENT)
                    if hasattr(subImage, "getMask"):
                        data = numpy.ma.masked_where(numpy.logical_or(subImage.getMask().getArray() & maskVal,
                                                                      subPoly.getArray() < 1.0),
                                                     subImage.getImage().getArray())
                    else:
                        data = numpy.ma.masked_where(subPoly.getArray() < 1.0, subImage.getArray())
                    data = data.compressed()
                    num = len(data)
                    value = data.mean() if num > 0 else numpy.nan
                else:
                    result = afwMath.makeStatistics(subImage, statistic | afwMath.NPOINT, stats)
                    value = result.getValue(statistic)
                    num = result.getValue(afwMath.NPOINT)

                assert(i < self._values.getWidth() and j < self._values.getHeight())

                self._values.set(i, j, self._values.get(i, j) + value*num)
                self._numbers.set(i, j, self._numbers.get(i, j) + num)

    @classmethod
    def fromImage(cls, config, image, box=None, polygon=None):
        if hasattr(image, "getMaskedImage"):
            image = image.getMaskedImage()
        if box is None:
            box = image.getBBox(afwImage.PARENT)
        self = cls(config, box=box, polygon=polygon)
        self.addImage(image)
        return self

    def merge(self, bg):
        if (self.config.xSize, self.config.ySize) != (bg.config.xSize, bg.config.ySize):
            raise RuntimeError("Size mismatch: %s vs %s" % ((self.config.xSize, self.config.ySize),
                                                            (bg.config.xSize, bg.config.ySize)))
        if self.box != bg.box:
            raise RuntimeError("Box mismatch: %s vs %s" % (self.box, bg.box))
        if self.polygon != bg.polygon:
            raise RuntimeError("Polygon mismatch: %s vs %s" % (self.polygon, bg.polygon))
        self._values += bg._values
        self._numbers += bg._numbers
        return self

    def __iadd__(self, bg):
        return self.merge(bg)

    def getStatsImage(self):
        values = self._values.clone()
        values /= self._numbers
        thresh = self.config.minFrac*self.config.xSize*self.config.ySize
        values.getArray()[:] = numpy.where(self._numbers.getArray() < thresh, numpy.nan, values.getArray())
        return self.approximateStatsImage(values)

    def approximateStatsImage(self, image):
        # XXX Doesn't handle NANs in the image
        if self.config.approximateStyle == "NONE":
            return image
        assert(self.config.approximateStyle == "CHEBYSHEV")
        approx = afwMath.ApproximateControl(afwMath.ApproximateControl.CHEBYSHEV,
                                            self.config.xApproximate, self.config.yApproximate)
        mi = afwImage.makeMaskedImage(image)
        return afwMath.makeApproximate(self._xCenter, self._yCenter, mi, self.box, approx).getImage()

    def getImage(self, *args):
        # Reproduce API for BackgroundMI.getImageF
        box = None
        method = "AKIMA_SPLINE"
        args = list(args)
        if args and isinstance(args[0], afwGeom.Box2I):
            box = args.pop(0)
        if args and isinstance(args[0], basestring):
            method = args.pop(0)
        # XXX leftovers?

        if box is None:
            box = self.box
        method = afwMath.stringToInterpStyle(method)

        values = self.getStatsImage().getArray()
        mask = numpy.isnan(values)
        values = numpy.ma.masked_where(mask, values.astype(float))

        targetImage = afwImage.ImageF(box)

        if self.polygon is not None:
            if not self.polygon.overlaps(afwGeom.Box2D(box)):
                targetImage.set(0.0)
                return targetImage
            polyBox = afwGeom.Box2I(self.polygon.getBBox())
            polyBox.clip(box)
            subImage = afwImage.ImageF(targetImage, polyBox, afwImage.PARENT)
            target = subImage.getArray()
            box = polyBox
        else:
            target = targetImage.getArray()

        xMin, yMin = box.getMin()
        xMax, yMax = box.getMax()
        width, height = box.getDimensions()

        def doInterpolate(xSample, ySample, xInterp, method, maxExtrapolate):
            if len(xSample) == 0:
                return numpy.ones_like(xInterp)*numpy.nan
            try:
                interpolated = afwMath.makeInterpolate(xSample, ySample, method).interpolate(xInterp)
            except:
                if method == afwMath.Interpolate.CONSTANT:
                    # We've already tried the most basic interpolation and it failed
                    return numpy.ones_like(xInterp)*numpy.nan
                interpolated = doInterpolate(xSample, ySample, xInterp,
                                             afwMath.lookupMaxInterpStyle(len(xSample)), maxExtrapolate)
            if self.config.extrapolationMode == "NONE":
                interpolated = numpy.where(numpy.logical_or(xInterp < xSample.min(),
                                                            xInterp > xSample.max()),
                                           self.config.extrapolateValue, interpolated)
            if self.config.extrapolationMode == "CONSTANT":
                maskSample = numpy.ma.masked_where(numpy.isnan(ySample), ySample).compressed()
                lowValue = maskSample[0]
                highValue = maskSample[-1]
                interpolated = numpy.where(xInterp < xSample.min() - maxExtrapolate, lowValue, interpolated)
                interpolated = numpy.where(xInterp > xSample.max() - maxExtrapolate, highValue, interpolated)
            if self.config.extrapolationMode == "LIMITED":
                interpolated = numpy.where(numpy.logical_or(xInterp < xSample.min() - maxExtrapolate,
                                                            xInterp > xSample.max() + maxExtrapolate),
                                           self.config.extrapolateValue, interpolated)
            return interpolated

        xMaxExtrapolate = self.config.maxExtrapolateFrac*self.config.xSize
        yMaxExtrapolate = self.config.maxExtrapolateFrac*self.config.ySize

        # Interpolation in x
        temp = numpy.zeros((width, len(self._yCenter)))
        xPosFull = numpy.arange(xMin, xMax + 1, dtype=float)
        for y, yCen in enumerate(self._yCenter):
            # XXX opt: limit yCen based on reach of interpolator
            xPos = numpy.ma.masked_where(mask[y,:], self._xCenter).compressed()
            temp[:,y] = doInterpolate(xPos, values[y,:].compressed(), xPosFull, method, xMaxExtrapolate)

        tempMask = numpy.isnan(temp)
        temp = numpy.ma.masked_where(tempMask, temp)

        # Interpolation in y
        yPosFull = numpy.arange(yMin, yMax + 1, dtype=float)
        for x in range(xMin, xMax + 1):
            dx = x - xMin
            yPos = numpy.ma.masked_where(tempMask[dx,:], self._yCenter).compressed()
            target[:,dx] = doInterpolate(yPos, temp[dx,:].compressed(), yPosFull, method, yMaxExtrapolate)
        # Throw away anything outside the polygon
        if self.polygon is not None:
            polyImage = self.polygon.createImage(box)
            polyArray = polyImage.getArray()
            target[:] = numpy.where(polyArray >= 1.0, target, 0.0)
            target[:] = numpy.where(numpy.logical_and(polyArray > 0.0, polyArray < 1.0), numpy.nan, target)

        return targetImage

    getImageF = getImage


class PolygonBackground(object):
    def __init__(self, config, box, polygons=[], bgList=None, badPolygons=[]):
        self.config = config
        self.box = box
        if bgList is None:
            self.goodPolygons, self.badPolygons = self.threshPolygons(polygons)
            bgList = [Background(config, box, poly) for poly in self.goodPolygons]
        else:
            self.goodPolygons = polygons
            self.badPolygons = badPolygons
        assert(len(bgList) == len(self.goodPolygons))
        self._bgList = bgList

    def __reduce__(self):
        return self.__class__, (self.config, self.box, self.goodPolygons, self._bgList, self.badPolygons)

    def addImage(self, image):
        for bg in self._bgList:
            bg.addImage(image)

    def threshPolygons(self, polygons):
        box = afwGeom.Box2D(self.box)
        minArea = self.config.minFrac*self.config.xSize*self.config.ySize
        minPerimeter = 2*self.config.xSize + 2*self.config.ySize
        goodPolygonList = []
        badPolygonList = []
        for poly in polygons:
            if not poly.overlaps(box):
                #print "Polygon %s not in box %s" % (poly, box)
                continue
            if poly.calculateArea() < minArea:
                #print "Polygon %s area is %f" % (poly, poly.calculateArea())
                badPolygonList.append(poly)
                continue
            if poly.calculatePerimeter() < minPerimeter:
                #print "Polygon %s perimeter is %f" % (poly, poly.calculatePerimeter())
                badPolygonList.append(poly)
                continue
            goodPolygonList.append(poly)
        return goodPolygonList, badPolygonList

    @classmethod
    def fromImage(cls, config, image, box=None, polygons=[]):
        if hasattr(image, "getMaskedImage"):
            image = image.getMaskedImage()
        if box is None:
            box = image.getBBox(afwImage.PARENT)
        self = cls(config, box, polygons)
        self.addImage(image)
        return self

    def merge(self, bg):
        if not self.config.compare(bg.config):
            raise RuntimeError("Config mismatch: %s vs %s" % (self.config, bg.config))
        if self.box != bg.box:
            raise RuntimeError("Box mismatch: %s vs %s" % (self.box, bg.box))

        for i, iPoly in enumerate(bg.goodPolygons):
            found = False
            for j, jPoly in enumerate(self.goodPolygons):
                if iPoly == jPoly:
                    self._bgList[j].merge(bg._bgList[i])
                    found = True
                    break
            if not found:
                self.goodPolygons.append(iPoly)
                self._bgList.append(bg._bgList[i])

        newBadPolygons = []
        if iPoly in bg.badPolygons:
            found = False
            for jPoly in self.goodPolygons:
                if iPoly == jPoly:
                    found = True
                    break
            for jPoly in self.badPolygons:
                if iPoly == jPoly:
                    found = True
                    break
            if not found:
                newBadPolygons.append(iPoly)
        self.badPolygons += newBadPolygons

        return self

    def __iadd__(self, bg):
        return self.merge(bg)

    def getStatsImage(self):
        image = None
        for bg in self._bgList:
            bgImage = bg.getStatsImage()
            if image is None:
                image = bgImage
            else:
                image += bgImage
        return image

    def getImage(self, *args):
        # Reproduce API for BackgroundMI.getImageF
        box = None
        method = "AKIMA_SPLINE"
        args = list(args)
        if args and isinstance(args[0], afwGeom.Box2I):
            box = args.pop(0)
        if args and isinstance(args[0], basestring):
            method = args.pop(0)
        # XXX leftovers?

        if box is None:
            box = self.box

        targetImage = afwImage.ImageF(box)
        targetImage.set(0.0)
        for bg in self._bgList:
            targetImage += bg.getImage(box, method)

        targetArray = targetImage.getArray()
        for poly in self.badPolygons:
            polyImage = poly.createImage(box)
            targetArray[:] = numpy.where(polyImage.getArray() > 0, numpy.nan, targetArray)

        return targetImage

    getImageF = getImage


if __name__ == "__main__":
    # make a ramping image (spline should be exact for linear increasing image
    nx = 512
    ny = 1024
    x0, y0 = 9876, 54321
    box = afwGeom.Box2I(afwGeom.Point2I(x0, y0), afwGeom.Extent2I(nx, ny))
    ramp = afwImage.ImageF(box)
    dzdx, dzdy, z0 = 0.5, 0.7, 10000.0
    dz2dx2, dz2dy2 = 0.01, -0.01

    for x in range(nx):
        for y in range(ny):
            ramp.set(x, y, dzdx*x + dzdy*y + z0 + dz2dx2*x**2 + dz2dy2*y**2)

    config = BackgroundConfig()
    config.xSize = 64
    config.ySize = 64

    bg = Background.fromImage(config, ramp)

    import lsst.afw.display.ds9 as ds9
    frame = 1
    ds9.mtv(ramp, frame=frame); frame += 1
    for interp in ("CONSTANT", "LINEAR", "NATURAL_SPLINE", "AKIMA_SPLINE"):
        diff = bg.getImageF(interp)
        ds9.mtv(diff, frame=frame); frame += 1
        diff -= ramp
        ds9.mtv(diff, frame=frame); frame += 1
        print interp, diff.getArray().sum(), diff.getArray().mean(), diff.getArray().std()

    stats = afwMath.StatisticsControl()
    stats.setNanSafe(True)
    bgCtrl = afwMath.BackgroundControl("AKIMA_SPLINE", ramp.getWidth()//config.xSize+1,
                                       ramp.getHeight()//config.ySize+1,
                                       "REDUCE_INTERP_ORDER", stats, "MEAN")
    afwBg = afwMath.cast_BackgroundMI(afwMath.makeBackground(ramp, bgCtrl))

    for interp in ("CONSTANT", "LINEAR", "NATURAL_SPLINE", "AKIMA_SPLINE"):
        diff = afwBg.getImageF(interp)
        ds9.mtv(diff, frame=frame); frame += 1
        diff -= ramp
        ds9.mtv(diff, frame=frame); frame += 1
        print "afw-" + interp, diff.getArray().sum(), diff.getArray().mean(), diff.getArray().std()

    ds9.mtv(bg.getStatsImage(), frame=frame); frame += 1
    ds9.mtv(afwBg.getStatsImage(), frame=frame); frame += 1

    def assertSame(image1, image2):
        if not numpy.all(image1.getArray() == image2.getArray()):
            raise AssertionError("Images don't match: %s %s" % (image1.getArray(), image2.getArray()))

    # Test pickle
    import pickle
    new = pickle.loads(pickle.dumps(bg))
    assertSame(bg.getImage(), new.getImage())

    # Check creation of sub-image
    box = afwGeom.Box2I(afwGeom.Point2I(123 + x0, 45 + y0), afwGeom.Extent2I(45, 123))
    bgImage = bg.getImageF()
    bgSubImage = afwImage.ImageF(bgImage, box, afwImage.PARENT)
    testImage = bg.getImageF(box, "AKIMA_SPLINE")
    assertSame(testImage, bgSubImage)
