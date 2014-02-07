#!/usr/bin/env python
'''
Solves Constrained Toy Problem Storing Optimization History.

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
import pdb

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Optimization
from pyOpt import SLSQP


# =============================================================================
# 
# =============================================================================
def objfunc(x):
	
	f = x[0]**2 + x[1]**2
	g = [0.0]*2
	g[0] = 3 - x[0]
	g[1] = 2 - x[1]
	
	fail = 0
	
	return f,g,fail
	

# =============================================================================
# 
# =============================================================================

# Instanciate Optimization Problem 
opt_prob = Optimization('TOY Constrained Problem',objfunc)
opt_prob.addVar('x1','c',value=1.0,lower=0.0,upper=10.0)
opt_prob.addVar('x2','c',value=1.0,lower=0.0,upper=10.0)
opt_prob.addObj('f')
opt_prob.addCon('g1','i')
opt_prob.addCon('g2','i')
print opt_prob

# Instanciate Optimizer (SLSQP) & Solve Problem Storing History
slsqp = SLSQP()
slsqp.setOption('IFILE','slsqp1.out')
slsqp(opt_prob,store_hst=True)
print opt_prob.solution(0)

# Solve Problem Using Stored History (Warm Start)
slsqp.setOption('IFILE','slsqp2.out')
slsqp(opt_prob, store_hst=True, hot_start='slsqp1')
print opt_prob.solution(1)
