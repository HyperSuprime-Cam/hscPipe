#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.processExposure import ProcessExposureTask
if __name__ == "__main__":
    submitPbs(ProcessExposureTask, "Reduce exposures using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processExposure.py'))

