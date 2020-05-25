#!/usr/bin/env python

import os
import sys

from .pyOpt_history import History
from .pyOpt_parameter import Parameter
from .pyOpt_variable import Variable
from .pyOpt_gradient import Gradient
from .pyOpt_constraint import Constraint
from .pyOpt_objective import Objective
from .pyOpt_optimization import Optimization
from .pyOpt_optimizer import Optimizer

__all__ = ['History','Parameter','Variable','Gradient','Constraint','Objective','Optimization','Optimizer']

dir = os.path.dirname(os.path.realpath(__file__))
for f in os.listdir(dir):
    if f.startswith('py') and os.path.isdir(os.path.join(dir,f)) and f not in ("pyIPOPT", ):
        try:
            exec('from .%s import %s' %(f, f[2:]))
            __all__.extend(sys.modules['pyOpt.'+f].__all__)
        except Exception as e:
            continue
