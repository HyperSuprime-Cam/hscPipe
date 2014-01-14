#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2014 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import numpy
import unittest
import lsst.utils.tests as utilsTests

import lsst.afw.coord as afwCoord
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage

from hsc.pipe.tasks.polygon import Polygon


def circle(radius, num, x0=0.0, y0=0.0):
    """Generate points on a circle

    @param radius: radius of circle
    @param num: number of points
    @param x0,y0: Offset in x,y
    @return x,y coordinates as numpy array
    """
    theta = numpy.linspace(0, 2*numpy.pi, num=num, endpoint=False)
    x = radius*numpy.cos(theta) + x0
    y = radius*numpy.sin(theta) + y0
    return numpy.array([x, y]).transpose()

class PolygonTest(unittest.TestCase):
    def setUp(self):
        self.x0 = 0.0
        self.y0 = 0.0

    def polygon(self, num, radius=1.0, x0=None, y0=None):
        """Generate a polygon

        @param num: Number of points
        @param radius: Radius of polygon
        @param x0,y0: Offset of center
        @return polygon
        """
        if x0 is None:
            x0 = self.x0
        if y0 is None:
            y0 = self.y0
        points = circle(radius, num, x0=x0, y0=y0)
        return Polygon([afwGeom.Point2D(x,y) for x,y in points])

    def square(self, size=1.0, x0=0, y0=0):
        """Generate a square

        @param size: Half-length of the sides
        @param x0,y0: Offset of center
        """
        return Polygon([afwGeom.Point2D(size*x + x0, size*y + y0)
                             for x, y in ((-1, -1), (-1, 1), (1, 1), (1, -1))])

    def testCenter(self):
        """Test Polygon.calculateCenter"""
        for num in range(3, 30):
            poly = self.polygon(num)
            center = poly.calculateCenter()
            self.assertAlmostEqual(center.getX(), self.x0)
            self.assertAlmostEqual(center.getY(), self.y0)

    def testContains(self):
        """Test Polygon.contains"""
        radius = 1.0
        for num in range(3, 30):
            poly = self.polygon(num, radius=radius)
            self.assertTrue(poly.contains(afwGeom.Point2D(self.x0, self.y0)))
            self.assertFalse(poly.contains(afwGeom.Point2D(self.x0 + radius, self.y0 + radius)))

    def testIntersects(self):
        """Test Polygon.intersects"""
        radius = 1.0
        for num in range(3, 30):
            poly1 = self.polygon(num, radius=radius)
            poly2 = self.polygon(num, radius=radius, x0=0.5*radius, y0=0.5*radius)
            poly3 = self.polygon(num, radius=2*radius)
            self.assertTrue(poly1.intersects(poly2))
            self.assertTrue(poly2.intersects(poly1))
            self.assertFalse(poly1.intersects(poly3))
            self.assertFalse(poly3.intersects(poly1))

    def testIntersectSquare(self):
        """Test Polygon.intersects with vertical edges (infinite slope)"""
        poly1 = self.square()
        poly2 = self.square(x0=1, y0=1)
        poly3 = self.square(x0=3, y0=3)
        self.assertTrue(poly1.intersects(poly2))
        self.assertTrue(poly2.intersects(poly1))
        self.assertFalse(poly1.intersects(poly3))
        self.assertFalse(poly3.intersects(poly1))

    def testOverlaps(self):
        """Test Polygon.overlaps"""
        return
        radius = 1.0
        for num in range(3, 30):
            poly1 = self.polygon(num, radius=radius)
            poly2 = self.polygon(num, radius=radius, x0=radius, y0=radius)
            poly3 = self.polygon(num, radius=2*radius)
            poly4 = self.polygon(num, radius=radius, x0=3*radius, y0=3*radius)
            self.assertTrue(poly1.overlaps(poly2))
            self.assertTrue(poly2.overlaps(poly1))
            self.assertTrue(poly1.overlaps(poly3))
            self.assertTrue(poly3.overlaps(poly1))
            self.assertFalse(poly1.overlaps(poly4))
            self.assertFalse(poly4.overlaps(poly1))

def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(PolygonTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
