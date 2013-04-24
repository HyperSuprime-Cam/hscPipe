#!/usr/bin/env python

import re
import os
import glob
import tempfile
import multiprocessing
from argparse import ArgumentParser

import pyfits

FILTERS = "grizy"

def subsetSchlafly(inName, outName):
    inFile = pyfits.open(inName)
    inData = inFile[1].data

    schema = pyfits.ColDefs([pyfits.Column(name="id", format="K"),
                             pyfits.Column(name="ra", format="D"),
                             pyfits.Column(name="dec", format="D")] +
                            [pyfits.Column(name=name, format="E") for name in FILTERS] +
                            [pyfits.Column(name=name + "_err", format="E") for name in FILTERS] +
                            [pyfits.Column(name=name + "_stdev", format="E") for name in FILTERS]
                            )

    outHdu = pyfits.new_table(schema, nrows=len(inData))
    outData = outHdu.data

    outData.ident = inData.obj_id
    outData.ra = inData.ra
    outData.dec = inData.dec
    for i, f in enumerate(FILTERS):
        # Some of the below (e.g., "mean") are functions in the pyfits.FITS_rec class,
        # so we need to access them differently than just grabbing an attribute.
        mean = outData.field(f)
        err = outData.field(f + "_err")
        stdev = outData.field(f + "_stdev")

        mean[:] = inData.field("mean")[:,i]
        err[:] = inData.field("err")[:,i]
        stdev[:] = inData.field("stdev")[:,i]

    outHdu.writeto(outName, clobber=True)
    print "Wrote %s" % outName
    inFile.close()

def system(command):
    print command
    os.system(command)


def hpsplit(inList, outroot, nside=16):
    args = "-r ra -d dec -n %d" % nside
    system("hpsplit -o " + outroot + " " + args + " " + " ".join(inList))

def generateIndexes(inName, outName, index, healpix=None, nside=None):
    args = "-S r -L 20 -E -M -j 0.4 -n 100"
    if healpix is not None:
        args += " -H %d" % healpix
    if nside is not None:
        args += " -s %d" % nside
    system("build-index -i " + inName + " -o " + outName + "_0.fits -I " + str(index) + "0 -P 0 " + args)
    system("build-index -1 " + outName + "_0.fits -o " + outName + "_1.fits -I " + str(index) + "1 -P 1 " + args)
    system("build-index -1 " + outName + "_0.fits -o " + outName + "_2.fits -I " + str(index) + "2 -P 2 " + args)
    system("build-index -1 " + outName + "_0.fits -o " + outName + "_3.fits -I " + str(index) + "3 -P 3 " + args)
    system("build-index -1 " + outName + "_0.fits -o " + outName + "_4.fits -I " + str(index) + "4 -P 4 " + args)


class _Caller(object):
    def __init__(self, func):
        self.func = func
    def __call__(self, args):
        return self.func(*args[0], **args[1])
class MapFunc(object):
    """Class to wrap the map() function with optional threading"""
    def __init__(self, threads, func, maxtasksperchild=1):
        self.func = func
        self.args = []
        if threads > 1:
            self.pool = multiprocessing.Pool(threads, maxtasksperchild=maxtasksperchild)
            self.map = self.pool.map
        else:
            self.map = map
    def add(self, *args, **kwargs):
        self.args.append((args, kwargs))
    def run(self):
        caller = _Caller(self.func)
        return self.map(caller, self.args)

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", nargs="*", help="Input files")
    parser.add_argument("-j", dest="threads", type=int, default=0, help="Number of threads")
    parser.add_argument("-o", "--output", required=True, help="Output root name")
    parser.add_argument("-s", "--nside", default=16, help="HEALPix nside (power of 2)")
    args = parser.parse_args()

    ps1List = []
    mapSubset = MapFunc(args.threads, subsetSchlafly)
    for i, inName in enumerate(args.input):
        ps1Name = "%s_in_%d.fits" % (args.output, i)
        ps1List.append(ps1Name)
        if not os.path.exists(ps1Name):
            mapSubset.add(inName, ps1Name)
    mapSubset.run()
    del mapSubset

    hpsplit(ps1List, "%s_hp_%%i.fits" % args.output, nside=args.nside)

    mapIndexes = MapFunc(args.threads, generateIndexes)
    for inName in glob.glob("%s_hp_*.fits" % args.output):
        m = re.search(r"%s_hp_(\d+)\.fits" % args.output, inName)
        assert m, "Unable to match filename"
        healpix = int(m.group(1))
        andName = "%s_and_%d" % (args.output, healpix)
        mapIndexes.add(inName, andName, healpix, nside=args.nside, healpix=healpix)
    mapIndexes.run()
