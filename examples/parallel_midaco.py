#!/usr/bin/env python
'''
Solves MIDACO's Toy Problem Using MIDACO's Parallelization

    min 	(x1-1.0)**2 + (x2-2.0)**2 + (x3-3.0)**2 + (x4-4.0)**2 + 1.23456789
    s.t.:	x1 - 1.0 = 0
            1.333333333 - x2 <= 0
            2.666666666 - x3 <= 0
            1 <= xi <= 4,  i = 1,...,4
        
    x* = [, ]
    f* = 
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, time

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
from pyOpt import MIDACO


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f =  (x[0]-1.0)**2 + (x[1]-2.0)**2 + (x[2]-3.0)**2 + (x[3]-4.0)**2 + 1.23456789
    
    g = [0.0]*3
    g[0] = x[0] - 1.0
    g[1] = 1.333333333 - x[1]
    g[2] = 2.666666666 - x[2]
    
    time.sleep(0.005)
    
    fail = 0
    return f,g, fail
    

# =============================================================================
# 
# =============================================================================

# Instantiate Optimization Problem 
opt_prob = Optimization('MIDACO Toy Problem',objfunc)
opt_prob.addVar('x1','c',lower=1.0,upper=4.0,value=1.0)
opt_prob.addVar('x2','c',lower=1.0,upper=4.0,value=1.0)
opt_prob.addVar('x3','c',lower=1.0,upper=4.0,value=1.0)
opt_prob.addVar('x4','c',lower=1.0,upper=4.0,value=1.0)
opt_prob.addObj('f')
opt_prob.addCon('g1','e')
opt_prob.addCon('g2','i')
opt_prob.addCon('g3','i')

# Solve Problem (No-Parallelization)
midaco_none = MIDACO()
midaco_none.setOption('IPRINT',-1)
midaco_none.setOption('MAXEVAL',50000)
midaco_none(opt_prob)
if myrank == 0:
    print opt_prob.solution(0)
#end

# Solve Problem (SPM-Parallelization)
midaco_spm = MIDACO(pll_type='SPM')
midaco_spm.setOption('IPRINT',-1)
midaco_none.setOption('MAXEVAL',50000)
midaco_spm(opt_prob)
print opt_prob.solution(1)
