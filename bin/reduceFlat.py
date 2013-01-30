#!/usr/bin/env python
import os
from hsc.pipe.base.pbs import submitPbs
from hsc.pipe.tasks.detrends import FlatTask
if __name__ == "__main__":
    submitPbs(FlatTask, "Construct flat using PBS",
              os.path.join(os.environ['HSCPIPE_DIR'], 'bin', 'processFlat.py'))
