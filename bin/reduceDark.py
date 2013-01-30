#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.detrends import DarkTask
if __name__ == "__main__":
    submitPbs(DarkTask, "Construct dark using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processDark.py'))
