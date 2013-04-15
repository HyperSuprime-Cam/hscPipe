import lsst.afw.image as afwImage
import lsst.afw.geom as afwGeom

def getTractPatchList(dataRef, skyMap):
    md = dataRef.get("calexp_md", immediate=True)
    wcs = afwImage.makeWcs(md)
    width, height = md.get("NAXIS1"), md.get("NAXIS2")
    pointList = [(0,0), (width, 0), (width, height), (0, height)]
    coordList = [wcs.pixelToSky(afwGeom.Point2D(x, y)) for x, y in pointList]
    tractPatchList = skyMap.findTractPatchList(coordList)
    # Convert from tract/patch objects to identifiers
    tractDict = dict((tractInfo.getId(), [patchInfo.getIndex() for patchInfo in patchList])
                     for tractInfo, patchList in tractPatchList)
    return tractDict

def invert(dataIdList, tractDictList):
    tractData = {}
    for dataId, tractDict in zip(dataIdList, tractDictList):
        for tractId, patchList in tractDict.items():
            if tractId not in tractData:
                tractData[tractId] = {}
            patchDict = tractData[tractId]
            for patchId in patchList:
                if patchId not in patchDict:
                    patchDict[patchId] = []
                patchDict[patchId].append(dataId)
    return tractData


def printDataRef(dataRef, tractDict, doPatches=True):
    print "CalExp %s: (#=%d)" % (dataRef.dataId, len(tractDict))
    for tract, patchList in tractDict.items():
        print "==> Tract %d (#=%d)" % (tract, len(patchList))
        if doPatches:
            print "====> Patches: %s" % patchList

def printTracts(tractData, tract=None, doDataId=False):
    if tract is not None:
        _printTract(tract, tractData[tract], doDataId=doDataId)
        return
    for tract in tractData:
        _printTract(tract, tractData[tract], doDataId=doDataId)

def _printTract(tract, patchDict, doDataId=False):
    print "Tract %d" % tract
    for patch, dataIdList in patchDict.items():
        print "==> Patch %s (#=%d)" % (patch, len(dataIdList))
        if doDataId:
            print "====> DataIds: %s" % dataIdList

def printPatch(tract, patch, tractData, doDataId=False):
    print "Tract %d, patch %s: (%s)" % (tract, patch, len(tractData[tract][patch]))
    if doDataId:
        print "==> DataIds: %s" % tractData[tract][patch]


__all__ = [getTractPatchList, invert, printDataRef, printTracts, printPatch]
