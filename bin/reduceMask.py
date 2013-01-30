#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.detrends import MaskTask
if __name__ == "__main__":
    submitPbs(MaskTask, "Construct mask using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processMask.py'))
