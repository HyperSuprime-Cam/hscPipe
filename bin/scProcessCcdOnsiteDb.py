#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import argparse, os, sys
from lsst.pipe.base import HscArgumentParser

from hsc.pipe.tasks.processCcd import SuprimeCamProcessCcdTask as TaskClass

# FH added for QA output
# import hsc.onsite.qa.fitsthumb as QaFitsthumb
# import hsc.onsite.qa.measSeeingQa as QaSeeing
#import hsc.onsite.qa.onsiteDbUtils as onsiteDbUtils
#import onsiteDbUtils as onsiteDbUtils
#import hsc.hscDb.frame_regist_CorrSuprime  as registCorrSup

if __name__ == "__main__":
    parser = HscArgumentParser(conflict_handler='resolve') # old style
    parser.add_argument('--dumpconfig', action="store_true", help="Dump the configuration to stdout and exit")

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        if "--doraise" in sys.argv:
            raise
        print >> sys.stderr, e
        sys.exit(1)

    if namespace.dumpconfig:
        namespace.config._save(sys.stdout)
        sys.exit(0)


    ## === update onsite Db status
    # !!! it is better to db update for frame_analysis_start just before execution of this script
    #     but, to get 'id' on time when this analysis process is invoked, I'm temporarily 
    #     doing this here. 
    def getDataId(namespace):  
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

    id, visit, ccd =  getDataId(namespace)
    onsiteDbUtils.updateStatusFrameAnalysisStart(id)
    
    ## === create task objects and run tasks 
    task = TaskClass(config=namespace.config)
    if False: ### debugging
        print '************ Here, config start **************'
        print namespace.config
        print '************ Here, config end **************'

    
    for sensorRef in namespace.dataRefList:
        
        if namespace.doRaise:
            task.run(sensorRef)
        else:
            if False:
                try:
                    task.run(sensorRef)
                except Exception, e:
                    task.log.log(task.log.FATAL, "Failed on dataId=%s: %s" % (sensorRef.dataId, e))
            else:
                print '* * '*40
                print sensorRef
                print '* * '*40
                task.run(sensorRef)

    ## === update onsite Db status
    onsiteDbUtils.updateStatusFrameAnalysisEnd(id)

    ## === register CORR data QA values
    filenameList = []
    for sensorRef in namespace.dataRefList:
        filenameList.append(sensorRef.get('calexp_filename')[0])
    ## for onsite Qa, only 1 file is input, so n of element should be 1
    filename = filenameList[0]
    print '**** CORR filename:', filename
    #registCorrSup.registCorrFrameMetaInfo(filename)
    onsiteDbUtils.registFrameQaMetaInfo(id, filename)
   
