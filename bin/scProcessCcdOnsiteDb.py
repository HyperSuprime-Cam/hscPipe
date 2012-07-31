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
from hsc.pipe.tasks import plotSetup # needed to enable non-gui matplotlib when DISPLAY is not set

from hsc.pipe.tasks.onsiteDb import SubaruProcessCcdOnsiteDbTask

# Note from Jim: the code that was previously here has been moved
# into the task class in python/hsc/pipe/tasks/onsiteDb.py, so we
# don't have to duplicate the parse-and-run code in pipe_base.
# This means the database is notified and filled each time a new
# CCD is processed, rather than once for all CCDs.  But a comment
# in the original code said that this script is only run on one
# CCD at a time, so I think it doesn't matter.

SubaruProcessCcdOnsiteDbTask.parseAndRun()