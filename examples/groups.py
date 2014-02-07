#!/usr/bin/env python
'''
Solves Constrained Toy Problem Using Variable Groups.

    min 	x1^2 + x2^2
    s.t.:	3 - x1 <= 0
            2 - x2 <= 0
            -10 <= x1 <= 10
            -10 <= x2 <= 10 
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Optimization
from pyOpt import SLSQP


# =============================================================================
# 
# =============================================================================
def objfunc(xn):
    
    x0 = xn['x']
    x1 = xn['z']
    
    f = x0**2 + x1**2
    g = [0.0]*2
    g[0] = 3 - x0
    g[1] = 2 - x1
    
    fail = 0
    
    return f,g,fail
    

# =============================================================================
# 
# ============================================================================= 

# Instantiate Optimization Problem
opt_prob = Optimization('TOY Constrained Problem',objfunc,use_groups=True)
opt_prob.addVarGroup('a',2,'c',value=1.0, lower=0.0, upper=10)
opt_prob.delVarGroup('a')
opt_prob.addVar('x','c',value=1.0, lower=0.0, upper=10)
opt_prob.addVarGroup('y',2,'c',value=1.0, lower=0.0, upper=10)
opt_prob.delVarGroup('y')
opt_prob.addVarGroup('z',1,'c',value=1.0, lower=0.0, upper=10)
opt_prob.addVarGroup('b',5,'c',value=3.0, lower=0.0, upper=10)
opt_prob.delVarGroup('b')
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
opt_prob.addCon('g2','i')
print opt_prob

# Instantiate Optimizer (SLSQP) & Solve Problem
slsqp = SLSQP()
slsqp(opt_prob)
print opt_prob.solution(0)

