#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.detrends import BiasTask
if __name__ == "__main__":
    submitPbs(BiasTask, "Construct bias using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processBias.py'))
