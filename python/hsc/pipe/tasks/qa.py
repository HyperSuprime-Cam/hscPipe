from lsst.pex.config import Config, ConfigurableField, Field, DictField
from lsst.pipe.base import Task
from . import measSeeingQa
from . import sizeMagnitudeMitakaStarSelector as starSel
import os
import numpy
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as cameraGeom


class QaConfig(Config):
    ##seeing = ConfigurableField(target=measSeeingQa.MeasureSeeingTask, doc="Measure seeing")
    seeing = ConfigurableField(target=measSeeingQa.MeasureSeeingMitakaTask, doc="Measure seeing ver. Mitaka")
    useIcsources = Field(dtype=bool, default=False, doc="Use icsources(calib.sources) rather than final sources") 
    doWriteMeta = Field(dtype=bool, default=False, doc="Write Qa metadata for CCD as FITS")
    magzeroExpected = DictField(keytype=str, itemtype=float, 
                                default={'g': 29.0, 'r': 29.0, 'i': 28.6, 'z': 27.7, 'y': 27.4},
                                doc="Expected magnzero (e/sec) in the case of ideal weather condition: for calcTransp")
    flatScale = DictField(keytype=str, itemtype=float, 
                                default={'g': 1.0, 'r': 1.0, 'i': 1.0, 'z': 1.0, 'y': 1.0},
                                doc="Count levels of the FOV center of flat field in each band: for calcTransp")

class QaTask(Task):
    ConfigClass = QaConfig
    _DefaultName = 'Qa'

    def __init__(self, *args, **kwargs):
        super(QaTask, self).__init__(*args, **kwargs)
        self.makeSubtask("seeing", metadata=self.metadata)

        ## initialization of QA values
        # common
        self.metadata.set('ANAID', -9999)
        self.metadata.set('REGISTID', -9999)
        self.metadata.set('VISIT', -9999)
        self.metadata.set('CCD_REGISTRY', -9999)
        self.metadata.set('ANAPATH', 'NOT_SET')
        self.metadata.set('RERUN', 'NOT_SET')
        self.metadata.set('FLAG_AUTO', -9999)
        self.metadata.set('FLAG_USR', -9999)
        self.metadata.set('FLAG_TAG', -9999)
        # seeing/starsel
        self.metadata.set('SEEING_MODE', -9999.0)
        self.metadata.set('ELL_MED', -9999.0)
        self.metadata.set('ELL_PA_MED', -9999.0)
        self.metadata.set("magLim", 99.0)
        self.metadata.set("TRANSP", -9999.0)
        self.metadata.set("fwhmRough", -9999.0)
        self.metadata.set("fwhmRobust", -9999.0)
        self.metadata.set("ellRobust", -9999.0)
        self.metadata.set("ellPaRobust", -9999.0)
        self.metadata.set("medianFwhmPsfSeq", -9999.0)
        self.metadata.set("sigmaFwhmPsfSeq", -9999.0)
        self.metadata.set("minFwhmPsfSeq", -9999.0)
        self.metadata.set("maxFwhmPsfSeq", -9999.0)
        self.metadata.set("numFwhmPsfLikeRobust", -9999)

        # astrometric match
        self.metadata.set("REFCAT", 'NOT_SET')

    def run(self, dataRef, exposure, sources):
        # initial set of exposure metadata
        #metadata = exposure.getMetadata()
        #metadata.combine(self.metadata)

        self.seeing.run(dataRef, sources, exposure)

        self.calcTransp(dataRef, exposure)

        # = info for qa process and result
        # config
        self.setRefcatName()

        # flags
        # this part should be done by calculating merit functions somewhere else in a polite manner.
        self.metadata.set('FLAG_AUTO', 0)
        self.metadata.set('FLAG_USR', 0)
        self.metadata.set('FLAG_TAG', 1)

        # = info for data management

        # registry information
        try:
            anaPath = dataRef.getButler().mapper.root
            dataId = dataRef.dataId
            if dataId.has_key('anaId') and dataId['anaId'] >= 0:
                anaId = dataId['anaId']
            else:
                anaId = -9999
            if dataId.has_key('registId') and dataId['registId'] >= 0:
                registId = dataId['registId']
            else:
                registId = -9999
            visit = dataId['visit']
            ccd = dataId['ccd']
            self.metadata.set('ANAID', anaId)
            self.metadata.set('REGISTID', registId)
            self.metadata.set('VISIT', visit)
            self.metadata.set('CCD_REGISTRY', ccd)
            self.metadata.set('ANAPATH', anaPath)

            self.log.info("registry info dataRef.dataId: %s (visit, ccd, registId, anaId): %s" % (str(dataId), str((visit, ccd, registId, anaId))))

        except Exception, e:
            self.log.warn("Could not get registry info (visit, ccd): %s" % str(e))

        try:
            rerunName = self.getRerunName(dataRef)
            self.metadata.set('RERUN', rerunName)
        except:
            self.log.warn("Could not determine rerun from output directory")

        # merging local Qa metadata into exposure metadata
        metadata = exposure.getMetadata()
        for key in self.metadata.names():
            metadata.set(key, self.metadata.get(key))

        if self.config.doWriteMeta:
            # preparing exposure object which holds exposureQA metadata
            expQaMeta = afwImage.ExposureI(0,0)
            expQaMeta.setCalib(exposure.getCalib())
            expQaMeta.setMetadata(metadata)
            self.log.info("writing an QA metadata for ccd %d at %s" % (ccd, dataRef.get('ccdMetadata_filename')[0]))
            try:
                dataRef.put(expQaMeta, 'ccdMetadata')
            except Exception, e:
                self.log.warn("Could not write an QA metadata for ccd: %s" % str(e))


    def getRerunName(self, sensorRef):
        """ derive rerun name from the file path, assuming the parent directory 'xxx/rerun' precedes the rerun name """
        # rerun: assuming directory tree has the fixed form where 'rerun/' is just followed by '$rerun_name/'
        corrPath = sensorRef.get('calexp_filename')[0]
        return corrPath[corrPath.find('rerun'):].split('/')[1]

    def setRefcatName(self):
        """ derive astrometry.net index name from env.variable, and set it to metadata. """
        try:
            catalogName = os.path.basename(os.environ.get('ASTROMETRY_NET_DATA_DIR').rstrip('/'))
            if catalogName == 'None': # for safety
                catalogName = 'NOT_SET'
        except:
            catalogName = 'NOT_SET'
        self.metadata.set('REFCAT', catalogName)


    def calcTransp(self, dataRef, exposure):
        """
        Procedure of estimating transparency based on magzero, setting values in related keywords
        """

        # Check if photocal succeeded
        try:
            fluxmag0 = exposure.getCalib().getFluxMag0()[0]
            if fluxmag0 <= 0:
                raise
        except Exception, e:
            self.log.warn("calcTransp: valid FLUXMAG0 is not present. Set to -9999.0: %s" % str(e))
            exposure.getCalib().setFluxMag0(-9999.0, -9999.0)
            self.metadata.set('MAGZERO', 99.0)
            transp = -9999.0
            self.log.info('calcTransp: transparency(estimated) is set to %f' % transp)
            return

        self.log.info("calcTransp: FLUXMAG0: %f" % fluxmag0)

        # magzero is available, then proceed.
        filter = exposure.getFilter().getName()
        self.log.info('calcTransp: filter = %s' % filter)

        metadata = exposure.getMetadata()
        try:
            # gain from amp info
            det = exposure.getDetector()
            ccd = cameraGeom.cast_Ccd(det)
            gainList = []
            for a in ccd:
                gain = a.getElectronicParams().getGain()
                gainList.append(gain)
            gain = numpy.mean(gainList)
            self.log.info('calcTransp: gainList: %s - meanGain = %f' % (str(gainList), gain))
            if False:
                gain = metadata.get('GAIN')
                if gain == 0.0: # HSC
                    gain1 = metadata.get('T_GAIN1')
                    gain2 = metadata.get('T_GAIN2')
                    gain3 = metadata.get('T_GAIN3')
                    gain4 = metadata.get('T_GAIN4')
                    gain = (gain1 + gain2 + gain3 + gain4)/4.0

            # scaling magzero's to (mag/e/sec)
            flatScale =  self.config.flatScale[filter]
            magzeroAduMeasured = metadata.get('MAGZERO') # mag/adu/sec
            magzeroMeasured = magzeroAduMeasured + 2.5*numpy.log10(gain * flatScale) # mag/e/sec
            magzeroExpected = self.config.magzeroExpected[filter] # mag/e/sec

            self.log.info('calcTransp: fluxscale used in flatfield = %f' % flatScale)
            self.log.info('calcTransp: magzero(measured) = %f (mag/e/s)' % magzeroMeasured)
            self.log.info('calcTransp: magzero(expected) = %f (mag/e/s)' % magzeroExpected)
            transp = (10.0**(0.4 * magzeroMeasured)) / (10.0**(0.4 * magzeroExpected))
            self.log.info('calcTransp: debug:transparency(estimated) = %f' % transp)
        except Exception, e:
            self.log.warn("Could not get gain or magzero: %s" % str(e))
            transp = -9999.0
            self.log.info('calcTransp: transparency(estimated) = %f' % transp)

        self.metadata.set('TRANSP', transp)

