import numpy as np
import lsst.pex.config
from lsst.pipe.tasks.processCoadd import ProcessCoaddTask
import lsst.daf.persistence as dafPersistence
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom

class PropagatePSFFlagsToCoaddConfig(ProcessCoaddTask.ConfigClass):
    tol = lsst.pex.config.Field(dtype=float, default=0.2, doc="Tolerance of matching in arcsec")
    propagatedFlags = lsst.pex.config.ListField(dtype=str, default=["calib.psf.used", "calib.psf.candidate"], doc="List of flags to be propagated")

class PropagatePSFFlagsToCoaddTask(ProcessCoaddTask):
    ConfigClass = PropagatePSFFlagsToCoaddConfig

    def run(self, dataRef):
        src = self.propagate(dataRef)
        self.write(dataRef, src)

    def propagate(self, dataRef):
        propagatedFlags = self.config.propagatedFlags

        # get source catalog
        src = dataRef.get("deepCoadd_src", immediate = True)
 
        # make new source catalog with new columns for PSF flags
        mapper = afwTable.SchemaMapper(src.schema)
        mapper.addMinimalSchema(src.schema)
        newSchema = mapper.getOutputSchema()
        keys = dict()
        for flag in propagatedFlags:
            for cond in ["any", "all"]:
                key = flag + "." + cond
                newSchema.addField(afwTable.Field["Flag"](key, "Source has %s = True in %s of matched objects in individual exposures" % (flag, cond)))
                keys[key] = newSchema.find(key).key
        newSrc = afwTable.SourceCatalog(newSchema)
        newSrc.extend(src, mapper=mapper)
        # Copy slots. After Jim pushes a new branch, we can remove this.
        for name in ("Centroid", "Shape", "ApFlux", "ModelFlux", "PsfFlux", "InstFlux"):
            meas = getattr(src.table, "get" + name + "Key")()
            err = getattr(src.table, "get" + name + "ErrKey")()
            flag = getattr(src.table, "get" + name + "FlagKey")()
            getattr(newSrc.table, "define" + name)(meas, err, flag)

        # get visit and CCD information
        coadd = dataRef.get("deepCoadd", immediate = True)
        ccdInputs = coadd.getInfo().getCoaddInputs().ccds
        visit = ccdInputs.get("visit")
        ccd = ccdInputs.get("ccd")
        
        # prepare butler for reding sources in indivisual visits
        butler = dataRef.getButler()
        srcWithFlag = dict()
        for flag in propagatedFlags:
            srcWithFlag[flag] = list()
        tol = self.config.tol/3600./180*np.pi
        closest = True

        # get flags from each visit
        for v in set(visit):
            srcWithFlagTmp = dict()
            for flag in propagatedFlags:
                srcWithFlagTmp[flag] = np.zeros(len(src), dtype = bool)
            for c in ccd[visit == v]:
                srcCcd = butler.get("src", dataId = {"visit": int(v), "ccd":int(c)}, immediate = True)
                for flag in propagatedFlags:
                    srcCcdWithFlag = srcCcd[srcCcd.get(flag)]
                    match = afwTable.matchRaDec(src[src.get('deblend.nchild') == 0], srcCcdWithFlag, afwGeom.Angle(tol), closest)
                    for m in match:
                        index = (np.where(src.get("id") == m.first.get("id")))[0][0]
                        srcWithFlagTmp[flag][index] = True
            for flag in propagatedFlags:
                srcWithFlag[flag].append(srcWithFlagTmp[flag])

        # collapse flags
        for flag in propagatedFlags:
            for cond in ["any", "all"]:
                collapsedSrcWithFlag = getattr(np, cond)(srcWithFlag[flag], axis = 0)
                for index in np.where(collapsedSrcWithFlag)[0]:
                    newSrc[index].setFlag(keys[flag+"."+cond], True)

        return newSrc

    def write(self, dataRef, src):
        dataRef.put(src, "deepCoadd_src")

