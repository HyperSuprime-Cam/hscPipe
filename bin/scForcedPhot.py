#!/usr/bin/env python
from hsc.pipe.tasks.deprecated import deprecated
deprecated("scForcedPhot.py", "hscForcedPhot.py")
from hsc.pipe.tasks.forcedPhot import SubaruForcedPhotTask
SubaruForcedPhotTask.parseAndRun()
