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
from lsst.pipe.base import ArgumentParser

from hsc.pipe.tasks.processCcd import SuprimeCamProcessCcdTask as TaskClass

# FH added for QA output
# import hsc.onsite.qa.fitsthumb as QaFitsthumb
# import hsc.onsite.qa.measSeeingQa as QaSeeing
#import hsc.onsite.qa.onsiteDbUtils as onsiteDbUtils
#import onsiteDbUtils as onsiteDbUtils
#import hsc.hscDb.frame_regist_CorrSuprime  as registCorrSup

class OutputAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if namespace.rerun:
            raise argparse.ArgumentTypeError("Please specify --output or --rerun, but not both")

        namespace.outPath = values

class RerunAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        """We can't just parse the arguments and reset namespace.outPath as the Mapper's been
        instantiated before we get a chance"""

        if namespace.outPath:
            raise argparse.ArgumentTypeError("Please specify --output or --rerun, but not both")

        envar = "SUPRIME_DATA_DIR"
        if os.environ.has_key(envar):
            namespace.rerun = values
            if namespace.camera == 'suprimecam':
                cameraKey = 'SUPA'
            elif namespace.camera == 'hsc':
                cameraKey = 'HSC'
            elif namespace.camera == 'hscsim':
                cameraKey = 'HSC'
            namespace.outPath = os.path.join(os.environ[envar], cameraKey, "rerun", namespace.rerun)
            if not os.path.exists(namespace.outPath):
                try:
                    os.makedirs(namespace.outPath) # should be in butler
                except OSError, e:
                    print "*** making output directories failed once."
                    if not e.errno == errno.EEXIST:
                        # if directory does not exist, something wrong occured
                        #raise RuntimeError, "Failed to create output directories: %s" % namespace.outPath
                        raise
        else:
            raise argparse.ArgumentTypeError("You must define $%s to use --rerun XXX" % envar)


if __name__ == "__main__":

    ## === parsing arguments 
    if sys.argv[1] == 'suprimecam':
        try:
            parser = ArgumentParser(name="suprimecam", conflict_handler='resolve') # new style
        except TypeError:
            parser = ArgumentParser(conflict_handler='resolve') # old style
    elif sys.argv[1] == 'hscsim':
        try:
            parser = ArgumentParser(name="hscsim", conflict_handler='resolve') # new style
        except TypeError:
            parser = ArgumentParser(conflict_handler='resolve') # old style
    else:
        try:
            parser = ArgumentParser(name="suprimecam", conflict_handler='resolve') # new style
        except TypeError:
            parser = ArgumentParser(conflict_handler='resolve') # old style

    parser.add_argument('--output', type=str, dest="outPath", default=None, help="output root directory",
                        action=OutputAction)
    parser.add_argument('--rerun', type=str, default=None, help='Desired rerun (overrides --output)',
                        action=RerunAction)

    try:
        namespace = parser.parse_args(config=TaskClass.ConfigClass())
    except Exception, e:
        print >> sys.stderr, e
        sys.exit(1)

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
   
