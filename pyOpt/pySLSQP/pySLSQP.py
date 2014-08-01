#!/usr/bin/env python
'''
pySLSQP - A Python pyOpt interface to SLSQP.

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.2   $Date: 21/06/2010 21:00$


Tested on:
---------
Linux with g77
Linux with gfortran
Linux with pathf95
Win32 with g77
Mac with g95

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2008)
	v. 1.1	- History support (PJ,RP, 2010)
	v. 1.2	- Gradient Class Support (PJ,RP, 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- store correct nevals from nfunc and ngrad
'''

# =============================================================================
# SLSQP Library
# =============================================================================
try:
	import slsqp
except ImportError:
	raise ImportError('SLSQP shared library failed to import')
#end

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import copy, time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Optimizer
from pyOpt import Gradient

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity
# =============================================================================
eps = 1.0	# define a value for machine precision
while ((eps/2.0 + 1.0) > 1.0):
	eps = eps/2.0
#end
eps = 2.0*eps
#eps = math.ldexp(1,-52)


# =============================================================================
# SLSQP Optimizer Class
# =============================================================================
class SLSQP(Optimizer):
	
	'''
	SLSQP Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		SLSQP Optimizer Class Initialization
		
		**Keyword arguments:**
		
		- pll_type -> STR: Parallel Implementation (None, 'POA'-Parallel Objective Analysis), *Default* = None
		
		Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
		'''
		
		#
		if (pll_type == None):
			self.poa = False
		elif (pll_type.upper() == 'POA'):
			self.poa = True
		else:
			raise ValueError("pll_type must be either None or 'POA'")
		#end
		
		#
		name = 'SLSQP'
		category = 'Local Optimizer'
		def_opts = {
		# SLSQP Options
		'ACC':[float,1e-6],			# Convergence Accurancy
		'MAXIT':[int,50], 			# Maximum Iterations
		'IPRINT':[int,1],			# Output Level (<0 - None, 0 - Screen, 1 - File)
		'IOUT':[int,6],     		# Output Unit Number
		'IFILE':[str,'SLSQP.out'],	# Output File Name
		}
		informs = {
		-1 : "Gradient evaluation required (g & a)",
		0 : "Optimization terminated successfully.",
		1 : "Function evaluation required (f & c)",
		2 : "More equality constraints than independent variables",
		3 : "More than 3*n iterations in LSQ subproblem",
		4 : "Inequality constraints incompatible",
		5 : "Singular matrix E in LSQ subproblem",
		6 : "Singular matrix C in LSQ subproblem",
		7 : "Rank-deficient equality constraint subproblem HFTI",
		8 : "Positive directional derivative for linesearch",
		9 : "Iteration limit exceeded",
		}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		
	def __solve__(self, opt_problem={}, sens_type='FD', store_sol=True, disp_opts=False, store_hst=False, hot_start=False, sens_mode='', sens_step={}, *args, **kwargs):
		
		'''
		Run Optimizer (Optimize Routine)
		
		**Keyword arguments:**
		
		- opt_problem -> INST: Optimization instance
		- sens_type -> STR/FUNC: Gradient type, *Default* = 'FD' 
		- store_sol -> BOOL: Store solution in Optimization class flag, *Default* = True 
		- disp_opts -> BOOL: Flag to display options in solution text, *Default* = False
		- store_hst -> BOOL/STR: Flag/filename to store optimization history, *Default* = False
		- hot_start -> BOOL/STR: Flag/filename to read optimization history, *Default* = False
		- sens_mode -> STR: Flag for parallel gradient calculation, *Default* = ''
		- sens_step -> FLOAT: Sensitivity setp size, *Default* = {} [corresponds to 1e-6 (FD), 1e-20(CS)]
		
		Additional arguments and keyword arguments are passed to the objective function call.
		
		Documentation last updated:  February. 2, 2011 - Peter W. Jansen
		'''
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pySLSQP - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pySLSQP: Parallel objective Function Analysis requires mpi4py'
			#end
			comm = MPI.COMM_WORLD
			nproc = comm.Get_size()
			if (mpi4py.__version__[0] == '0'):
				Bcast = comm.Bcast
			elif (mpi4py.__version__[0] == '1'):
				Bcast = comm.bcast
			#end
			self.pll = True
			self.myrank = comm.Get_rank()
		else:
			self.pll = False
			self.myrank = 0
		#end
		
		myrank = self.myrank
		
		#
		def_fname = self.options['IFILE'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		#
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)
		
		
		#======================================================================
		# SLSQP - Objective/Constraint Values Function
		#======================================================================
		def slfunc(m,me,la,n,f,g,x):
			
			# Variables Groups Handling
			if opt_problem.use_groups:
				xg = {}
				for group in group_ids.keys():
					if (group_ids[group][1]-group_ids[group][0] == 1):
						xg[group] = x[group_ids[group][0]]
					else:
						xg[group] = x[group_ids[group][0]:group_ids[group][1]]
					#end
				#end
				xn = xg
			else:
				xn = x
			#end
			
			# Flush Output Files
			self.flushFiles()
			
			# Evaluate User Function
			fail = 0
			ff = []
			gg = []
			if (myrank == 0):
				if self.h_start:
					[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						[ff,gg,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
					#end
				#end
			#end
			
			if self.pll:
				self.h_start = Bcast(self.h_start,root=0)
			#end
			if self.h_start and self.pll:
				[ff,gg,fail] = Bcast([ff,gg,fail],root=0)
			elif not self.h_start:	
				[ff,gg,fail] = opt_problem.obj_fun(xn, *args, **kwargs)
			#end
			
			# Store History
			if (myrank == 0):
				if self.sto_hst:
					log_file.write(x,'x')
					log_file.write(ff,'obj')
					log_file.write(gg,'con')
					log_file.write(fail,'fail')
				#end
			#end	
			
			# Objective Assigment
			if isinstance(ff,complex):
				f = ff.astype(float)
			else:
				f = ff
			#end
			
			# Constraints Assignment (negative gg as slsqp uses g(x) >= 0)
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(gg[i],complex):
					g[i] = -gg[i].astype(float)
				else:
					g[i] = -gg[i]
				#end
			#end
			
			return f,g
		
		
		#======================================================================
		# SLSQP - Objective/Constraint Gradients Function
		#======================================================================
		def slgrad(m,me,la,n,f,g,df,dg,x):
			
			if self.h_start:
				dff = []
				dgg = []
				if (myrank == 0):
					[vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						dff = vals['grad_obj'][0].reshape((len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())))
						dgg = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))	
					#end
				#end
				if self.pll:
					self.h_start = Bcast(self.h_start,root=0)
				#end
				if self.h_start and self.pll:
					[dff,dgg] = Bcast([dff,dgg],root=0)
				#end
			#end
			
			if not self.h_start:	
				#
				dff,dgg = gradient.getGrad(x, group_ids, [f], -g, *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')
			#end
			
			# Gradient Assigment
			for i in xrange(len(opt_problem._variables.keys())):
				df[i] = dff[0,i]
				for jj in xrange(len(opt_problem._constraints.keys())):
					dg[jj,i] = -dgg[jj,i]
				#end
			#end
			
			return df,dg
		
		
		
		# Variables Handling
		n = len(opt_problem._variables.keys())
		xl = []
		xu = []
		xx = []
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'c'):
				xl.append(opt_problem._variables[key].lower)
				xu.append(opt_problem._variables[key].upper)
				xx.append(opt_problem._variables[key].value)
			elif (opt_problem._variables[key].type == 'i'):
				raise IOError('SLSQP cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('SLSQP cannot handle discrete design variables')
			#end
		#end
		xl = numpy.array(xl)
		xu = numpy.array(xu)
		xx = numpy.array(xx)
		
		# Variables Groups Handling
		group_ids = {}
		if opt_problem.use_groups:
			k = 0
			for key in opt_problem._vargroups.keys():
				group_len = len(opt_problem._vargroups[key]['ids'])
				group_ids[opt_problem._vargroups[key]['name']] = [k,k+group_len]
				k += group_len
			#end
		#end
		
		# Constraints Handling
		m = len(opt_problem._constraints.keys())
		meq = 0
		#gg = []
		if m > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					meq += 1
				#end
				#gg.append(opt_problem._constraints[key].value)
			#end
		#end
		#gg = numpy.array(gg,numpy.float)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff,numpy.float)
		
		
		# Setup argument list values
		la = numpy.array([max(m,1)], numpy.int)
		gg = numpy.zeros([la], numpy.float)
		n1 = numpy.array([n+1], numpy.int)
		df = numpy.zeros([n+1], numpy.float)
		dg = numpy.zeros([la,n+1], numpy.float)
		acc = numpy.array([self.options['ACC'][1]], numpy.float)
		maxit = numpy.array([self.options['MAXIT'][1]], numpy.int)
		iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
		if (myrank != 0):
			iprint = -1
		else:
			iprint = self.options['IPRINT'][1]
		#end
		iout = numpy.array([self.options['IOUT'][1]], numpy.int)
		ifile = self.options['IFILE'][1]
		if (iprint >= 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		mode = numpy.array([0], numpy.int)
		mineq = m - meq + 2*(n+1)
		lsq = (n+1)*((n+1)+1) + meq*((n+1)+1) + mineq*((n+1)+1)
		lsi = ((n+1)-meq+1)*(mineq+2) + 2*mineq
		lsei = ((n+1)+mineq)*((n+1)-meq) + 2*meq + (n+1)
		slsqpb = (n+1)*(n/2) + 2*m + 3*n + 3*(n+1) + 1
		lwM = lsq + lsi + lsei + slsqpb + n + m
		lw = numpy.array([lwM], numpy.int)
		w = numpy.zeros([lw], numpy.float)
		ljwM = max(mineq,(n+1)-meq)
		ljw = numpy.array([ljwM], numpy.int)
		jw = numpy.zeros([ljw], numpy.intc)
		nfunc = numpy.array([0], numpy.int)
		ngrad = numpy.array([0], numpy.int)
		
		# Run SLSQP
		t0 = time.time()
		slsqp.slsqp(m,meq,la,n,xx,xl,xu,ff,gg,df,dg,acc,maxit,iprint,
			iout,ifile,mode,w,lw,jw,ljw,nfunc,ngrad,slfunc,slgrad)
		sol_time = time.time() - t0
		
		if (myrank == 0):
			if self.sto_hst:
				log_file.close()
				if tmp_file:
					hos_file.close()
					name = hos_file.filename
					os.remove(name+'.cue')
					os.remove(name+'.bin')
					os.rename(name+'_tmp.cue',name+'.cue')
					os.rename(name+'_tmp.bin',name+'.bin')
				#end
			#end
		#end
		
		if (iprint > 0):
			slsqp.closeunit(self.options['IOUT'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = mode[0]
		sol_inform['text'] = self.getInform(mode[0])
		
		if store_sol:
			
			sol_name = 'SLSQP Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = 0
			
			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = xx[i]
				i += 1
			#end
			
			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = ff[i]
				i += 1
			#end
			
			if m > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = -gg[i]
					i += 1
				#end
			else:
				sol_cons = {}
			#end
			
			sol_lambda = {}
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, Sensitivities=sens_type,
				myrank=myrank, arguments=args, **kwargs)
			
		#end
			
		return ff, xx, sol_inform
		
		
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getOption(self, name):
		
		'''
		Get Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getInform(self, infocode):
		
		'''
		Get Optimizer Result Information (Optimizer Specific Routine)
		
		Keyword arguments:
		-----------------
		id -> STRING: Option Name
		
		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		'''
		
		return self.informs[infocode]
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		'''
		
		# 
		iPrint = self.options['IPRINT'][1]
		if (iPrint >= 0):
			slsqp.pyflush(self.options['IOUT'][1])	
		#end
	


#==============================================================================
# SLSQP Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test SLSQP
	print 'Testing ...'
	slsqp = SLSQP()
	print slsqp
	
