#!/usr/bin/env python

from lsst.pipe.tasks.run import runTask
from hsc.pipe.tasks.processCcd import SuprimeCamProcessCcdTask

runTask(SuprimeCamProcessCcdTask)

