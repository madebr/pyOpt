#!/usr/bin/env python
'''
Solves Schittkowski's TP37 Constrained Problem.

    min 	-x1*x2*x3
    s.t.:	x1 + 2.*x2 + 2.*x3 - 72 <= 0
            - x1 - 2.*x2 - 2.*x3 <= 0
            0 <= xi <= 42,  i = 1,2,3
    
    f* = -3456 , x* = [24, 12, 12]
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time
import pdb

# =============================================================================
# Extension modules
# =============================================================================
#from pyOpt import *
from pyOpt import Optimization
from pyOpt import PSQP
from pyOpt import SLSQP
from pyOpt import CONMIN
from pyOpt import COBYLA
from pyOpt import SOLVOPT
from pyOpt import KSOPT
from pyOpt import NSGA2
from pyOpt import ALGENCAN
from pyOpt import FILTERSD


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f = -x[0]*x[1]*x[2]
    g = [0.0]*2
    g[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    g[1] = -x[0] - 2.*x[1] - 2.*x[2]
    
    fail = 0
    return f,g, fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem 
opt_prob = Optimization('TP37 Constrained Problem',objfunc)
opt_prob.addVar('x1','c',lower=0.0,upper=42.0,value=10.0)
opt_prob.addVar('x2','c',lower=0.0,upper=42.0,value=10.0)
opt_prob.addVar('x3','c',lower=0.0,upper=42.0,value=10.0)
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
opt_prob.addCon('g2','i')
print opt_prob

# Instantiate Optimizer (PSQP) & Solve Problem
psqp = PSQP()
psqp.setOption('IPRINT',0)
psqp(opt_prob,sens_type='FD')
print opt_prob.solution(0)

# Instantiate Optimizer (SLSQP) & Solve Problem
slsqp = SLSQP()
slsqp.setOption('IPRINT',-1)
slsqp(opt_prob,sens_type='FD')
print opt_prob.solution(1)

# Instantiate Optimizer (CONMIN) & Solve Problem
conmin = CONMIN()
conmin.setOption('IPRINT',0)
conmin(opt_prob,sens_type='CS')
print opt_prob.solution(2)

# Instantiate Optimizer (COBYLA) & Solve Problem
cobyla = COBYLA()
cobyla.setOption('IPRINT',0)
cobyla(opt_prob)
print opt_prob.solution(3)

# Instantiate Optimizer (SOLVOPT) & Solve Problem
solvopt = SOLVOPT()
solvopt.setOption('iprint',-1)
solvopt(opt_prob,sens_type='FD')
print opt_prob.solution(4)

# Instantiate Optimizer (KSOPT) & Solve Problem
ksopt = KSOPT()
ksopt.setOption('IPRINT',0)
ksopt(opt_prob,sens_type='FD')
print opt_prob.solution(5)

# Instantiate Optimizer (NSGA2) & Solve Problem
nsga2 = NSGA2()
nsga2.setOption('PrintOut',0)
nsga2(opt_prob)
print opt_prob.solution(6)

# Instantiate Optimizer (ALGENCAN) & Solve Problem
algencan = ALGENCAN()
algencan.setOption('iprint',0)
algencan(opt_prob)
print opt_prob.solution(7)

# Instantiate Optimizer (FILTERSD) & Solve Problem
filtersd = FILTERSD()
filtersd.setOption('iprint',0)
filtersd(opt_prob)
print opt_prob.solution(8)

