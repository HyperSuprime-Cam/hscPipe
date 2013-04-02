#!/usr/bin/env python

import os
import sys

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base.struct import Struct
import lsst.afw.table as afwTable

from .onsite import SubaruProcessCcdOnsiteTask
from .processCcd import SubaruProcessCcdConfig

class OnsiteDbTask(pipeBase.Task):
    ConfigClass = pexConfig.Config
    
    def start(self, registId):
        self.onsiteDbUtils.updateStatusFrameAnalysisStart(registId)

    def end(self, registId):
        self.onsiteDbUtils.updateStatusFrameAnalysisEnd(registId)

    def frameMetadata(self, registId, filename):
        self.onsiteDbUtils.registFrameQaMetaInfo(registId, filename)

    def getRegisterId(self, sensorRef):
        # XXX currently unused?
        dataId = sensorRef.dataId
        ccd = int(dataId['ccd'])
        visit = int(dataId['visit'])
        return self.onsiteDbUtils.getRegistryId(visit, ccd)

class HscOnsiteDbTask(OnsiteDbTask):
    def __init__(self, *args, **kwargs):
        super(HscOnsiteDbTask, self).__init__(*args, **kwargs)
        import hsc.onsite.onsiteDbUtilsSuprime as onsiteDbUtils
        self.onsiteDbUtils = onsiteDbUtils

class SuprimecamOnsiteDbTask(OnsiteDbTask):
    def __init__(self, *args, **kwargs):
        super(HscOnsiteDbTask, self).__init__(*args, **kwargs)
        import hsc.onsite.onsiteDbUtilsHsc as onsiteDbUtils
        self.onsiteDbUtils = onsiteDbUtils


class SubaruProcessCcdOnsiteDbConfig(SubaruProcessCcdConfig):
    onsiteDb = pexConfig.ConfigurableField(target=OnsiteDbTask, doc="Task for onsite database interaction")

class SubaruProcessCcdOnsiteDbTask(SubaruProcessCcdOnsiteTask):
    """Subclass of SubaruProcessCcdOnsiteTask that uses the database.
    """
    ConfigClass = SubaruProcessCcdOnsiteDbConfig
    _DefaultName = "processCcdOnsiteDb"

    def __init__(self, *args, **kwargs):
        super(SubaruProcessCcdOnsiteDbTask, self).__init__(*args, **kwargs)
        self.makeSubtask("onsiteDb")

    @classmethod
    def _makeArgumentParser(cls):
        argumentParser = super(SubaruProcessCcdOnsiteDbTask, cls)._makeArgumentParser()
        argumentParser.add_argument("--anaid", action="store", type=int, dest="anaId",
                                    default=0, help="analysis session Id in analysis table.")
        argumentParser.add_argument("--registid", action="store", type=int, dest="registId",
                                    default=0, help="primary key Id in registry raw table.")
        argumentParser.add_argument("--qa-logdir", action="store", type=str, dest="qaLogDir",
                                    default='.', help="directory where logging files in QaTasks or DbAccess are stored.")
        return argumentParser

    def run(self, sensorRef):
        # !!! it is better to db update for frame_analysis_start just before execution of this script
        #     but, to get 'id' on time when this analysis process is invoked, I'm temporarily 
        #     doing this here.

        # XXX self.parsedCmd is not available in an upgraded lsst.pipe.base.CmdLineTask; we will need to update
        anaId =  self.parsedCmd.anaId
        registId = self.parsedCmd.registId
        sensorRef.dataId['anaId'] = anaId
        sensorRef.dataId['registId'] = registId 
        self.log.info("anaid: %d registid: %d" % (sensorRef.dataId['anaId'], sensorRef.dataId['registId']))

        qaLogDir = self.parsedCmd.qaLogDir
        
        self.onsiteDb.start(registId)
        # Run main processing task and QA by calling base class
        result = SubaruProcessCcdOnsiteTask.run(self, sensorRef)
        ## === update onsite Db status
        self.onsiteDb.end(registId)
        ## === register CORR data QA values
        filename = sensorRef.get('calexp_filename')[0]
        self.onsiteDb.frameMetadata(registId, filename)
        return result
