# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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


import numpy as np

class PhotometricMagnitude(object):
    """ Convert between photometric fluxes and magnitudes 
    
    Simplify the conversion. Create an object with a zero point 
    magnitude and flux, and do the conversion with function calls,
    instead of messing around with logs and factors of -2.5
    """
    
    def __init__(self, zeroFlux=1, zeroMag=0):
        """Initialise the object with the flux corresponding to 
        a given magnitude. 
        
        Example:
        If a known 12th magnitude object has a flux of 15000 in an image, 
        set zeroFlux=15000 and zeroMag to 12
        """
        
        self.zeroFlux = float(zeroFlux)
        self.zeroMag =  float(zeroMag)
        

    def __str__(self):
        return "PhotometricMagnitude: magnitude %g at %g" \
            %(self.zeroMag, self.zeroFlux)

    def __call__(self, flux):
        return self.getMag(flux)
                    
    def getMag(self, flux):
        """Convert a flux to a magnitude
        
            flux can be a scalar, a list, a numpy array, or
            anything else that can be converted into a numpy array
        """
        
        #Convert the input type to a numpy array, do the tranformation
        #then convert back to the original type
        inputType = type(flux)
        output = self.zeroMag - 2.5*np.log10(np.array(flux)/self.zeroFlux)
        
        #Convert back to input type if necessary
        if inputType != type(np.array([0])):
            output = inputType(output)

        return output
        
    
    def getFlux(self, mag):
        """Convert a magnitude to a flux 
        
            flux can be a scalar, a list, a numpy array, or
            anything else that can be converted into a numpy array
        """

        #Convert the input type to a numpy array, do the tranformation
        #then convert back to the original type
        inputType = type(mag)
        inputMag = np.array(mag)
        
        output = -1*(inputMag - self.zeroMag)/2.5
        output = self.zeroFlux * np.power(10, output)

        return inputType(output)
    
