#!/usr/bin/env python
"""Solves Schittkowski's TP37 Constrained Problem.

min         -x1*x2*x3
s.t.:       x1 + 2.*x2 + 2.*x3 - 72 <= 0
        - x1 - 2.*x2 - 2.*x3 <= 0
        0 <= xi <= 42,  i = 1,2,3

f* = -3456 , x* = [24, 12, 12]
"""

from pyOpt import Optimization


def objfunc(x):
    f = -x[0] * x[1] * x[2]
    g = [
        x[0] + 2. * x[1] + 2. * x[2] - 72.0,
        -x[0] - 2. * x[1] - 2. * x[2],
    ]

    fail = 0
    return f, g, fail


def getlastsolution(prob):
    new_index = prob.firstavailableindex(prob.getSolSet())
    return prob.getSol(new_index - 1)


# Instantiate Optimization Problem
opt_prob = Optimization('TP37 Constrained Problem', objfunc)
opt_prob.addVar('x1', 'c', lower=0.0, upper=42.0, value=10.0)
opt_prob.addVar('x2', 'c', lower=0.0, upper=42.0, value=10.0)
opt_prob.addVar('x3', 'c', lower=0.0, upper=42.0, value=10.0)
opt_prob.addObj('f')
opt_prob.addCon('g1', 'i')
opt_prob.addCon('g2', 'i')
print(opt_prob)

# Instantiate Optimizer (PSQP) & Solve Problem
from pyOpt import PSQP

psqp = PSQP()
psqp.setOption('IPRINT', 0)
psqp(opt_prob, sens_type='FD')
print(getlastsolution(opt_prob))

# Instantiate Optimizer (SLSQP) & Solve Problem
from pyOpt import SLSQP

slsqp = SLSQP()
slsqp.setOption('IPRINT', -1)
slsqp(opt_prob, sens_type='FD')
print(getlastsolution(opt_prob))

# Instantiate Optimizer (CONMIN) & Solve Problem
from pyOpt import CONMIN

conmin = CONMIN()
conmin.setOption('IPRINT', 0)
conmin(opt_prob, sens_type='CS')
print(getlastsolution(opt_prob))

# Instantiate Optimizer (COBYLA) & Solve Problem
from pyOpt import COBYLA

cobyla = COBYLA()
cobyla.setOption('IPRINT', 0)
cobyla(opt_prob)
print(getlastsolution(opt_prob))

# Instantiate Optimizer (SOLVOPT) & Solve Problem
from pyOpt import SOLVOPT

solvopt = SOLVOPT()
solvopt.setOption('iprint', -1)
solvopt(opt_prob, sens_type='FD')
print(getlastsolution(opt_prob))

# Instantiate Optimizer (KSOPT) & Solve Problem
from pyOpt import KSOPT

ksopt = KSOPT()
ksopt.setOption('IPRINT', 0)
ksopt(opt_prob, sens_type='FD')
print(getlastsolution(opt_prob))

# Instantiate Optimizer (NSGA2) & Solve Problem
from pyOpt import NSGA2

nsga2 = NSGA2()
nsga2.setOption('PrintOut', 0)
nsga2(opt_prob)
print(getlastsolution(opt_prob))

# Instantiate Optimizer (ALGENCAN) & Solve Problem
from pyOpt import ALGENCAN

algencan = ALGENCAN()
algencan.setOption('iprint', 0)
algencan(opt_prob)
print(getlastsolution(opt_prob))
#
# # Instantiate Optimizer (FILTERSD) & Solve Problem
from pyOpt import FILTERSD

filtersd = FILTERSD()
filtersd.setOption('iprint', 0)
filtersd(opt_prob)
print(getlastsolution(opt_prob))
