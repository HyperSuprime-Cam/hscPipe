#!/usr/bin/env python

import os
import base64
import hashlib
import glob
import sqlite3
import argparse
import cPickle as pickle

from lsst.daf.persistence import Butler

def md5(thing):
    """Return a readable MD5 hash of the pickled thing"""
    return base64.b16encode(hashlib.md5(pickle.dumps(thing, -1)).digest())

def imageRepr(image):
    """Get the scientific representation of the image

    We only really care about the pixels.  Everything else is in the registry,
    which we are verifying separately.
    """
    impl = lambda img: img.getArray()
    if hasattr(image, "getMaskedImage"):  # Exposure
        image = image.getMaskedImage()
    if hasattr(image, "getMask"):  # MaskedImage
        return (impl(image.getImage()), impl(image.getMask()), impl(image.getVariance()))
    if hasattr(image, "getImage"):  # DecoratedImage
        return impl(image.getImage())
    return impl(image)  # Image


def getCalibs(registry, calibType):
    """Pull a list of dataIds for all calibs from the registry"""
    columns = {"validStart": str,
               "validEnd": str,
               "calibDate": str,
               "filter": str,
               "calibVersion": str,
               "ccd": int
               }
    query = "SELECT %s from %s" % (",".join(columns.keys()), calibType)
    cursor = registry.cursor()
    values = cursor.execute(query)
    return [{key: columns[key](value) for key, value in zip(columns.keys(), row)} for row in values]


def verifyCalibs(butler, registry):
    for calibType in ("bias", "dark", "flat", "fringe"):
        dataIdList = getCalibs(registry, calibType)
        print calibType, md5(dataIdList)
        for dataId in dataIdList:
            print dataId, md5(imageRepr(butler.get(calibType, dataId, visit=0, immediate=True)))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Get checksums for calibs")
    parser.add_argument("root", help="Data repo root directory (for butler)")
    parser.add_argument("calibRoot", help="Calibs root directory")
    args = parser.parse_args()

    butler = Butler(args.root, calibRoot=args.calibRoot)
    registry = sqlite3.connect(os.path.join(args.calibRoot, "calibRegistry.sqlite3"))
    verifyCalibs(butler, registry)
