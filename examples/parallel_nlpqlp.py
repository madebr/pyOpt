#!/usr/bin/env python
'''
Solves Schittkowski's TP37 Constrained Problem Using NLPQLP's Parallel Line Search 

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
from pyOpt import NLPQLP


# =============================================================================
# 
# =============================================================================
def objfunc(x):
    
    f = -x[0]*x[1]*x[2]
    g = [0.0]*2
    g[0] = x[0] + 2.*x[1] + 2.*x[2] - 72.0
    g[1] = -x[0] - 2.*x[1] - 2.*x[2]
    
    time.sleep(0.5)
    
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

# Solve Problem (No-Parallelization)
nlpqlp_none = NLPQLP()
nlpqlp_none.setOption('IPRINT',0)
nlpqlp_none(opt_prob)
if myrank == 0:
    print opt_prob.solution(0)
#end

# Solve Problem (Parallel Gradient)
nlpqlp_pgc = NLPQLP()
nlpqlp_pgc.setOption('IPRINT',0)
nlpqlp_pgc(opt_prob,sens_mode='pgc')
if myrank == 0:
    print opt_prob.solution(1)
#end

# Solve Problem (Parallel Line Search)
nlpqlp_spm = NLPQLP(pll_type='SPM')
nlpqlp_spm.setOption('IPRINT',0)
nlpqlp_spm(opt_prob)
print opt_prob.solution(2)
