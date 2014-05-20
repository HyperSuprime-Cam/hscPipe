
import hashlib
import zlib
import cPickle as pickle

import lsst.afw.image as afwImage

class CheckSum(object):

    sumFunctions = {
        "CRC32" : lambda obj: zlib.crc32(pickle.dumps(obj)),
        "MD5"   : lambda obj: hashlib.md5(pickle.dumps(obj)).hexdigest()
    }
    
    def __init__(self, sumType='MD5'):
        self.reset(sumType)
        
    def reset(self, sumType):
        if sumType in self.sumFunctions:
            self.sumType = sumType
        else:
            raise ValueError("Sum type must be in: "+",".join(self.sumFunctions.keys()))

    def __call__(self, obj):
        return self.sumFunctions[self.sumType](obj)

        
class ImageCheckSum(CheckSum):

    maskedImageTypes = (
        afwImage.MaskedImageF, afwImage.MaskedImageD,
        afwImage.DecoratedImageF, afwImage.DecoratedImageD,
        )
    exposureTypes = (
        afwImage.ExposureF, afwImage.ExposureD,
        )
    imageTypes = (
        afwImage.ImageF, afwImage.ImageI, afwImage.ImageD,
        )
    
    def __init__(self, stride=8, **kwargs):
        self.stride = stride
        super(ImageCheckSum, self).__init__(**kwargs)
        
    def __call__(self, obj):
        
        if isinstance(obj, self.exposureTypes):
            mimg = obj.getMaskedImage()
            img = mimg.getImage()
        elif isinstance(obj, self.maskedImageTypes):
            img = obj.getImage()
        elif isinstance(obj, self.imageTypes):
            img = obj
        else:
            raise RuntimeError("ImageCheckSum only works on afw Image,MaskedImage, and Exposure objects.")
            
        nimg = img.getArray()[::self.stride,::self.stride]
        sum = super(ImageCheckSum, self).__call__(nimg)
        
        return sum
