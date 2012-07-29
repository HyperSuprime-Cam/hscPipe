#!/usr/bin/env python

import os

import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from .onsite import SubaruProcessCcdOnsiteTask

class SubaruProcessCcdOnsiteDbTask(SubaruProcessCcdOnsiteTask):
    """Subclass of SubaruProcessCcdOnsiteTask that uses the database.
    """

    def run(self, sensorRef):
        # !!! it is better to db update for frame_analysis_start just before execution of this script
        #     but, to get 'id' on time when this analysis process is invoked, I'm temporarily 
        #     doing this here. 
        self.importDbUtils()
        id, visit, ccd =  self.getDataId()
        self.onsiteDbUtils.updateStatusFrameAnalysisStart(id)
        # Run main processing task and QA by calling base class
        result = OnsiteSubaruProcessCcdTask.run(self, sensorRef)
        ## === update onsite Db status
        self.onsiteDbUtils.updateStatusFrameAnalysisEnd(id)
        ## === register CORR data QA values
        filename = sensorRef.get('calexp_filename')[0]
        #registCorrSup.registCorrFrameMetaInfo(filename)
        self.onsiteDbUtils.registFrameQaMetaInfo(id, filename)
        return result

    def importDbUtils(self):
        namespace = self.parsedCmd
        try:
            if namespace.camera in ['suprimecam', 'sc', 'suprimecam-mit', 'mit']:
                import onsiteDbUtilsSuprime as onsiteDbUtils
            elif namespace.camera in ['hsc', 'hscsim']:
                import onsiteDbUtilsHsc as onsiteDbUtils
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
        if namespace.camera in ['suprimecam', 'sc', 'suprimecam-mit', 'mit']:
            id = int(visit)*10 + int(ccd)
        elif namespace.camera in ['hsc', 'hscsim']:
            #### !!! TBD how to assign visit for HSC data
            id = int(visit)*1000 + int(ccd)
            #id = int(visit)*100 + int(ccd)
        else:
            print >> sys.stderr, "Given instrument name is not valid: %s" % namespace.camera
            sys.exit(1)
        return id, visit, ccd
