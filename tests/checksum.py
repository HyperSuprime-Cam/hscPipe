#!/usr/bin/env python

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

import unittest
import lsst.utils.tests as utilsTests

import lsst.afw.image as afwImage
import lsst.afw.geom  as afwGeom

import hsc.pipe.tasks.checksum as checksum

class CheckSumTest(unittest.TestCase):

    def setUp(self):
        
        nx, ny = 256, 256

        self.known = {
            # values for an afwImage initialized with zeros
            'MD5_IMG'  : '97a008e029be34390d479b114e641a6a',
            'CRC32_IMG': -660391333,
            
            # values for an afwImage.Mask initialized with zeros
            'MD5_MSK'  : '17d2c6554aefb075d4c24e5202cc365f',
            'CRC32_MSK': 814355089,

            # values for the native list object
            'MD5'      : 'a32f230cdf14635312fdd6ba0f75c7f6',
            'CRC32'    : -536728847,
        }
        self.known['MD5_VAR']   = self.known['MD5_IMG']
        self.known['CRC32_VAR'] = self.known['CRC32_IMG']

        
        self.img   = afwImage.ImageF(nx, ny)
        self.dimg  = afwImage.DecoratedImageF(afwGeom.Extent2I(nx, ny))
        self.mimg  = afwImage.MaskedImageF(nx, ny)
        self.exp   = afwImage.ExposureF(nx, ny)
        self.list  = range(10)
        
    def tearDown(self):
        del self.img, self.mimg, self.exp, self.dimg


    def testReset(self):
        """Test that reset() method works"""
        
        sumObj = checksum.CheckSum()
        for hash in "MD5", "CRC32":
            sumObj.reset(hash)
            label = hash + "_IMG"
            self.assertEqual(sumObj(self.img)[label], self.known[label])
                         
    def testTypes(self):
        """Test that checksum works on image types"""
        
        for x in self.img, self.mimg, self.exp, self.dimg, self.list:
            for hash in "MD5", "CRC32":
                sumObj = checksum.CheckSum(hash)
                sums = sumObj(x)
                for k, v in sums.iteritems():
                    print k, v, self.known[k]
                    self.assertEqual(v, self.known[k])

                    
def suite():
    """Returns a suite containing all the test cases in this module."""

    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(CheckSumTest)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)

