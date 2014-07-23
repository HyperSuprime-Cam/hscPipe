
import hashlib
import zlib
import cPickle as pickle

import lsst.afw.image as afwImage

class CheckSum(object):

    exposureTypes = (
        afwImage.ExposureF, afwImage.ExposureD,
        )
    maskedImageTypes = (
        afwImage.MaskedImageF, afwImage.MaskedImageD,
        )
    decoratedImageTypes = (
        afwImage.DecoratedImageF, afwImage.DecoratedImageD,
        )
    imageTypes = (
        afwImage.ImageF, afwImage.ImageD, afwImage.ImageI,
        )
    
    protocol = 2
    sumFunctions = {
        "CRC32" : lambda obj: zlib.crc32(pickle.dumps(obj, CheckSum.protocol)),
        "MD5"   : lambda obj: hashlib.md5(pickle.dumps(obj, CheckSum.protocol)).hexdigest(),
    }
    
    def __init__(self, sumType='MD5'):
        self.reset(sumType)
        
    def reset(self, sumType):
        if sumType in self.sumFunctions:
            self.sumType = sumType
        else:
            raise ValueError("Sum type must be in: "+",".join(self.sumFunctions.keys()))

    def __call__(self, obj):

        func = self.sumFunctions[self.sumType]
        sum = {}

        # handle maskedImages
        if isinstance(obj, (self.exposureTypes, self.maskedImageTypes)):
            if isinstance(obj, self.exposureTypes):
                mimg = obj.getMaskedImage()
            else:
                mimg = obj
            sum[self.sumType + "_IMG"] = func(mimg.getImage())
            sum[self.sumType + "_MSK"] = func(mimg.getMask())
            sum[self.sumType + "_VAR"] = func(mimg.getVariance())
            
        # handle Images
        elif isinstance(obj, (self.decoratedImageTypes, self.imageTypes)):
            if isinstance(obj, self.decoratedImageTypes):
                img = obj.getImage()
            else:
                img = obj
            sum[self.sumType + "_IMG"] = func(img)
            
        # handle everything else
        else:
            sum[self.sumType] = func(obj)
            
        return sum

        
