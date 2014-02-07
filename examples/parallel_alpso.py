#!/usr/bin/env python
'''
Solves Runarsson's G08 Problem Using ALPSO Parallelization Options

    min 	-(sin(2*pi*x1)**3*sin(2*pi*x2))
    s.t.:	x1**2 - x2 + 1 <= 0
            1 - x1 + (x2-4)**2 <= 0
            0 <= xi <= 10,  i = 1,2
        
    x* = [1.2279713, 4.2453733]
    f* = -0.095825
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time, math

# =============================================================================
# External Python modules
# =============================================================================
try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
except:
    raise ImportError('mpi4py is required for parallelization')
#end

# =============================================================================
# Extension modules
# =============================================================================
#from pyOpt import *
from pyOpt import Optimization
from pyOpt import ALPSO


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f = -(((math.sin(2*math.pi*x[0])**3)*math.sin(2*math.pi*x[1]))/((x[0]**3)*(x[0]+x[1])))
    
    g = [0.0]*2
    g[0] = x[0]**2 - x[1] + 1
    g[1] = 1 - x[0] + (x[1]-4)**2
    
    time.sleep(0.01)
    
    fail = 0
    return f,g, fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem 
opt_prob = Optimization('G08 Global Constrained Problem',objfunc)
opt_prob.addVar('x1','c',lower=5.0,upper=1e-6,value=10.0)
opt_prob.addVar('x2','c',lower=5.0,upper=1e-6,value=10.0)
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
opt_prob.addCon('g2','i')

# Solve Problem (No-Parallelization)
alpso_none = ALPSO()
alpso_none.setOption('fileout',0)
alpso_none(opt_prob)
if myrank == 0:
    print opt_prob.solution(0)
#end

# Solve Problem (SPM-Parallelization)
alpso_spm = ALPSO(pll_type='SPM')
alpso_spm.setOption('fileout',0)
alpso_spm(opt_prob)
print opt_prob.solution(1)

# Solve Problem (DPM-Parallelization)
alpso_dpm = ALPSO(pll_type='DPM')
alpso_dpm.setOption('fileout',0)
alpso_dpm(opt_prob)
print opt_prob.solution(2)
