#!/usr/bin/env python

import os
import sys

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
from lsst.pipe.base.struct import Struct
import lsst.afw.table as afwTable

from .onsite import SubaruProcessCcdOnsiteTask

class SubaruProcessCcdOnsiteDbTask(SubaruProcessCcdOnsiteTask):
    """Subclass of SubaruProcessCcdOnsiteTask that uses the database.
    """

    @classmethod
    def parseAndRun(cls, args=None, config=None, log=None):
        """
        Copied from pipe_base.cmdLineTask.CmdLineTask,
        to add extra options for analysis operation (anaid, registid).
        """
        argumentParser = cls._makeArgumentParser()
        argumentParser.add_argument("--anaid", action="store", type=int, dest="anaId",
                                    default=0, help="analysis session Id in analysis table.")
        argumentParser.add_argument("--registid", action="store", type=int, dest="registId",
                                    default=0, help="primary key Id in registry raw table.")

        if config is None:
            config = cls.ConfigClass()
        parsedCmd = argumentParser.parse_args(config=config, args=args, log=log, overrides=cls.overrides)
        task = cls(name = cls._DefaultName, config = parsedCmd.config, log = parsedCmd.log)
        task.parsedCmd = parsedCmd
        task.runDataRefList(parsedCmd.dataRefList, doRaise=parsedCmd.doraise)
        return Struct(
            argumentParser = argumentParser,
            parsedCmd = parsedCmd,
            task = task,
            )

    def run(self, sensorRef):
        # !!! it is better to db update for frame_analysis_start just before execution of this script
        #     but, to get 'id' on time when this analysis process is invoked, I'm temporarily 
        #     doing this here.
        self.importDbUtils()
        if False:
            registid, visit, ccd =  self.getDataId()

        anaId =  self.parsedCmd.anaId
        registId = self.parsedCmd.registId
        sensorRef.dataId['anaId'] = anaId
        sensorRef.dataId['registId'] = registId 
        self.log.info("anaid: %d registid: %d" % (sensorRef.dataId['anaId'], sensorRef.dataId['registId']))

        self.onsiteDbUtils.updateStatusFrameAnalysisStart(registId)
        # Run main processing task and QA by calling base class
        result = SubaruProcessCcdOnsiteTask.run(self, sensorRef)
        ## === update onsite Db status
        self.onsiteDbUtils.updateStatusFrameAnalysisEnd(registId)
        ## === register CORR data QA values
        filename = sensorRef.get('calexp_filename')[0]
        #registCorrSup.registCorrFrameMetaInfo(filename)
        self.onsiteDbUtils.registFrameQaMetaInfo(registId, filename)
        return result

    def importDbUtils(self):
        namespace = self.parsedCmd
        try:
            if namespace.camera.lower() in ['suprimecam', 'sc', 'suprimecam-mit', 'mit']:
                import hsc.onsite.onsiteDbUtilsSuprime as onsiteDbUtils
            elif namespace.camera.lower() in ['hsc', 'hscsim']:
                import hsc.onsite.onsiteDbUtilsHsc as onsiteDbUtils
            else:
                print >> sys.stderr, "Given instrument name is not valid: %s" % namespace.camera
                sys.exit(1)
        except Exception, e:
            print >> sys.stderr, e
            sys.exit(1)
        self.onsiteDbUtils = onsiteDbUtils  # stuff module in self so other methods can use it

    def getDataId(self):
        namespace = self.parsedCmd
        dataId = (namespace.dataIdList)[0]
        ccd = int(dataId['ccd'])
        visit = int(dataId['visit'])

        # id is obtained from the one place (sup/hscpipe db) rather than being calculated in individual places
        if False:
            if namespace.camera.lower() in ['suprimecam', 'sc', 'suprimecam-mit', 'mit']:
                id = int(visit)*10 + int(ccd)
            elif namespace.camera.lower() in ['hsc', 'hscsim']:
                #### !!! TBD how to assign visit for HSC data
                #id = int(visit)*1000 + int(ccd)  
                id = int(visit)*100 + int(ccd)
            else:
                print >> sys.stderr, "Given instrument name is not valid: %s" % namespace.camera
                sys.exit(1)
        else:
            registId = self.onsiteDbUtils.getRegistryId(visit, ccd)

        return registId, visit, ccd



