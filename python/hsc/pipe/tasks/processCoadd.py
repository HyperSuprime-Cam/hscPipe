from lsst.pex.config import Field
from lsst.pipe.tasks.processCoadd import ProcessCoaddTask
import hsc.pipe.base.matches as hscMatches

class SubaruProcessCoaddConfig(ProcessCoaddTask.ConfigClass):
    doWriteDenormalizedMatches = Field(dtype=bool, default=True, doc="Write the denormalized match table?")

class SubaruProcessCoaddTask(ProcessCoaddTask):
    ConfigClass = SubaruProcessCoaddConfig

    def run(self, dataRef):
        results = super(SubaruProcessCoaddTask, self).run(dataRef)
        self.write(dataRef, results)
        return results

    def write(self, dataRef, struct):
        prefix = self.config.coaddName + "Coadd_"
        if self.config.doWriteDenormalizedMatches:
            dataRef.put(hscMatches.matchesToCatalog(struct.matches, struct.matchMeta), prefix + "srcMatchFull")

