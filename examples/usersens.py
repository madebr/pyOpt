#!/usr/bin/env python
'''
Solves Constrained Rosen-Suzuki Function

    - Different gradient calculation approaches are used
    
    f* = -6 , x* = [0, 1, 2, -1]
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Optimization
from pyOpt import CONMIN


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f = x[0]**2 - 5.*x[0] + x[1]**2 - 5.*x[1] + 2.*x[2]**2 - 21.*x[2] + x[3]**2 + 7.0*x[3] + 50.
    g = [0.0]*3
    g[0] = x[0]**2 + x[0] + x[1]**2 - x[1] + x[2]**2 + x[2] + x[3]**2 - x[3] - 8.0
    g[1] = x[0]**2 - x[0] + 2. * x[1]**2 + x[2]**2 + 2.*x[3]**2 - x[3] - 10.0
    g[2] = 2.*x[0]**2 + 2.*x[0] + x[1]**2 - x[1] + x[2]**2 - x[3] -5.0
    
    fail = 0
    
    return f,g,fail
    

def gradfunc(x,f,g):
    
    g_obj = [0.0]*4
    g_obj[0] = 2.*x[0] - 5
    g_obj[1] = 2.*x[1] - 5
    g_obj[2] = 4.*x[2] - 21
    g_obj[3] = 2.*x[3] + 7
    
    g_con = numpy.zeros([3,4])
    g_con[0][0] = 2.*x[0] + 1
    g_con[0][1] = 2.*x[1] - 1
    g_con[0][2] = 2.*x[2] + 1
    g_con[0][3] = 2.*x[3] - 1
    g_con[1][0] = 2.*x[0] - 1
    g_con[1][1] = 4.*x[1]
    g_con[1][2] = 2.*x[2]
    g_con[1][3] = 4.*x[3] - 1
    g_con[2][0] = 4.*x[0] + 2
    g_con[2][1] = 2.*x[1] - 1
    g_con[2][2] = 2.*x[2]
    g_con[2][3] = -1.
    
    fail = 0
    
    return g_obj,g_con,fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem
opt_prob = Optimization('Constrained Rosen-Suzuki',objfunc)
opt_prob.addVar('x1','c',value=1.5)
opt_prob.addVar('x2','c',value=1.5)
opt_prob.addVar('x3','c',value=1.5)
opt_prob.addVar('x4','c',value=1.5)
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
opt_prob.addCon('g2','i')
opt_prob.addCon('g3','i')
print opt_prob

# Instantiate Optimizer (CONMIN)
conmin = CONMIN()

# Solve Problem with Optimizer Using Finite Differences
conmin(opt_prob,sens_type='FD')
print opt_prob.solution(0)

# Solve Problem with Optimizer Using Complex Step
conmin(opt_prob,sens_type='CS')
print opt_prob.solution(1)

# Solve Problem with Optimizer Using User-Provided Sensitivities
conmin(opt_prob,sens_type=gradfunc)
print opt_prob.solution(2)
