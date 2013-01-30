#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.detrends import FringeTask
if __name__ == "__main__":
    submitPbs(FringeTask, "Construct fringe using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processFringe.py'))
