#!/usr/bin/env python
'''
Solves Constrained Rosenbrock's Passing Arguments into Objective Function.

    min 	A1*(X(2)-X(1)^2)^2 + (A2-x(1))^2
    s.t.:	X(1)^2 + X(2)^2 - A2 <= 0
            -1.0 <= xi <= 1.0,  i = 1,2
            
    arguments (passing into objfunc as kwargs)
        A1 = 100.0, A2 = 1.0    (kwarg passing as list of args)
        A3 = 1.0                (kwarg passing as single arg)
    
    f* = 0.0457 , x* = [0.7864, 0.6177]
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
from pyOpt import SLSQP


# =============================================================================
# 
# =============================================================================
def objfunc(x, **kwargs):
    
    a1 = kwargs['a12'][0]
    a2 = kwargs['a12'][1]
    a3 = kwargs['a3']
    
    f = a1*(x[1]-x[0]**2.)**2. + (a2-x[0])**2.
    g = [0.0]*2
    g[0] = x[0]**2. + x[1]**2.0 - a3
    
    fail = 0
    return f,g, fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem 
opt_prob = Optimization('Rosenbrock Constrained Problem',objfunc)
opt_prob.addVar('x1','c',lower=0.0,upper=1.0,value=0.5)
opt_prob.addVar('x2','c',lower=0.0,upper=1.0,value=0.5)
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
print opt_prob

# Arguments to pass into objfunc
a1 = 100.0
a2 = 1.0
a3 = 1.0

# Instantiate Optimizer (SLSQP) & Solve Problem
slsqp = SLSQP()
slsqp.setOption('IPRINT',-1)
slsqp(opt_prob,sens_type='FD',a12=[a1,a2],a3=a3)
print opt_prob.solution(0)
