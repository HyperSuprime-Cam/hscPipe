#!/usr/bin/env python

"""
This module is intended to "monkey patch" some LSST codes for HSC.

WARNING: Monkey patching is software skullduggery that should only be used as a
last resort.  It involves hacking modules at run time, by a different project
than which originally delivered the code being modified.  Because it is
deceiving the user, it should be considered as a dirty, evil hack that is only
occasionally warranted.  See http://en.wikipedia.org/wiki/Monkey_patch .

Motivation
==========

The LSST packages pex_config and pipe_base continued to evolve after we forked
for HSC.  Unfortunately, these packages provide base classes that are used in
functional code that we care about for HSC.  If we want to update that
functional code, we will get changes that rely on particular API changes in
pex_config and pipe_base that we are not (yet?) willing to introduce into the
HSC stack (because it would cause a period of instability as everything gets
updated).  We therefore need to redirect functionality directed at code in
pex_config and pipe_base that we do not yet support.

Implementation
==============

The largest and most visible change in the LSST packages is the introduction of
lsst.pex.config.ConfigurableField, and its use in lsst.pipe.base.Task.  The
idea behind this is to be able to make run-time changes to the task hierarchy,
to dynamically change which algorithms are used (e.g., selecting an algorithm
appropriate for a particular camera).  We will emulate the behaviour of
ConfigurableField using ConfigField and a registry to generate the correct
subtask.

NOTE: The ConfigurableField.retarget() functionality (to make run-time changes
to which Task subclass should be used) is NOT reproduced.  This means users
wanting to change the Task hierarchy will need to do so statically by
subclassing all the way up the hierarchy, overloading the __init__ of each
Task that is subclassed to use self.makeSubtask() to select the correct
subclass of its children tasks.

We also make a bare-bones definition of lsst.pipe.base.CmdLineTask (intended
to make scripts be simple two-liners: import the task and call its parseAndRun()
method).  Our definition is merely to prevent failure due to the symbol not
being defined.  We do not provide an implementation of parseAndRun() as it's
not (yet?) part of the HSC way.

Use
===

To activate this, merely "import" this module before any other module that may
use ConfigurableField or CmdLineTask.  A 'NOTE' is printed to advise the user
that underhanded practises are afoot.
"""

import lsst.pex.config
import lsst.pipe.base

if not 'ConfigurableField' in dir(lsst.pex.config):
    global REGISTRY # Mapping from Config class to Task
    REGISTRY = {} 

    class MP_ConfigurableField(lsst.pex.config.ConfigField):
        """A monkey patch for lsst.pex.config.ConfigurableField

        ConfigurableField is defined in LSST code, but not in HSC code.
        It is replaced with ConfigField and a registry.
        """
        def __init__(self, doc, target, ConfigClass=None):
            """Associate a target lsst.pipe.base.Task subclass with its appropriate
            configuration class so the Task.makeSubtask() method can create the right
            Task.
            """
            if ConfigClass is None and hasattr(target, 'ConfigClass'):
                ConfigClass = target.ConfigClass
            if ConfigClass is None:
                raise RuntimeError("Unable to determine configuration class for %s" % target.__name__)
            lsst.pex.config.ConfigField.__init__(self, doc=doc, dtype=ConfigClass)
            REGISTRY[ConfigClass] = target

        def retarget(self, *args, **kwargs):
            raise NotImplementedError("This replacement ConfigurableField cannot retarget; " +
                                      "you'll have to subclass all the way up.")

    class MP_Task(lsst.pipe.base.Task):
        """A monkey patch for lsst.pipe.base.Task

        Task.makeSubtask in LSST code uses ConfigurableField,
        which is not in HSC code.  We use a registry to achieve
        similar behaviour, but without the dynamic retargeting.
        """
        def makeSubtask(self, name, *args, **kwargs):
            """Create a subtask as a new instance self.<name>"""
            config = getattr(self.config, name, None)
            if config is not None and config.__class__ in REGISTRY:
                # Need some sleight-of-hand to fiddle with the API
                TaskClass = REGISTRY[config.__class__]
                subtask = TaskClass(config=config, name=name, parentTask=self, **kwargs)
                setattr(self, name, subtask)
                return
            # Otherwise, nothing to see here; move along
            super(MP_Task, self).makeSubtask(name, *args, **kwargs)

    # Patch them in
    print "NOTE: Monkey-patching lsst.pex.config.ConfigurableField"
    setattr(lsst.pex.config, 'ConfigurableField', MP_ConfigurableField)
    setattr(lsst.pipe.base, 'Task', MP_Task)


if not 'CmdLineTask' in dir(lsst.pipe.base):
    class MP_CmdLineTask(lsst.pipe.base.Task):
        @classmethod
        def parseAndRun(cls, *args, **kwargs):
            raise NotImplementedError("This replacement CmdLineTask cannot parseAndRun; call run yourself")

    print "NOTE: Monkey-patching lsst.pipe.base.CmdLineTask"
    setattr(lsst.pipe.base, 'CmdLineTask', MP_CmdLineTask)
