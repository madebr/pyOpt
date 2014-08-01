#!/usr/bin/env python
'''
pyCOBYLA - A Python pyOpt interface to COBYLA.

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.3   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with g77
Linux with pathf95
Win32 with g77

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2008)
	v. 1.1	- Integrate to pyOpt Framework (RP, 2009)
	v. 1.2	- History support (PJ,RP, 2010)
	v. 1.3	- Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# COBYLA Library
# =============================================================================
try:
	import cobyla
except:
	raise ImportError('COBYLA shared library failed to import')
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
# COBYLA Optimizer Class
# =============================================================================
class COBYLA(Optimizer):
	
	'''
	COBYLA Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		COBYLA Optimizer Class Initialization
		
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
		name = 'COBYLA'
		category = 'Local Optimizer'
		def_opts = {
		'RHOBEG':[float,0.5],		# Initial Variables Change
		'RHOEND':[float,1.0e-6],	# Convergence Accurancy
		'IPRINT':[int,2],			# Print Flag (0 - None, 1 - Final, 2,3 - Iteration)
		'MAXFUN':[int,3500],     	# Maximum Iterations
		'IOUT':[int,6],     		# Output Unit Number
		'IFILE':[str,'COBYLA.out'],	# Output File Name
		}
		informs = {
		0: 'Normal return',
		1: 'Max. number of function evaluations reach',
		2: 'Rounding errors are becoming damaging',
		}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		
	def __solve__(self, opt_problem={}, store_sol=True, disp_opts=False, store_hst=False, hot_start=False, *args, **kwargs):
		
		'''
		Run Optimizer (Optimize Routine)
		
		**Keyword arguments:**
		
		- opt_problem -> INST: Optimization instance
		- store_sol -> BOOL: Store solution in Optimization class flag, *Default* = True 
		- disp_opts -> BOOL: Flag to display options in solution text, *Default* = False
		- store_hst -> BOOL/STR: Flag/filename to store optimization history, *Default* = False
		- hot_start -> BOOL/STR: Flag/filename to read optimization history, *Default* = False
		
		Additional arguments and keyword arguments are passed to the objective function call.
		
		Documentation last updated:  February. 2, 2011 - Peter W. Jansen
		'''
		
		#
		if self.poa:
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyCOBYLA: Parallel objective Function Analysis requires mpi4py'
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
		
		
		#======================================================================
		# COBYLA - Objective/Constraint Values Function
		#======================================================================
		def calcfc(n,m,x,f,g):
			
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
			
			# Evaluate User Function (Real Valued)
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
			
			# Constraints Assigment (negative gg as cobyla uses g(x) >= 0)
			i = 0
			for j in xrange(len(opt_problem._constraints.keys())):
				if isinstance(gg[j],complex):
					g[i] = -gg[j].astype(float)
				else:
					g[i] = -gg[j]
				#end
				i += 1
			#end
			j = 0
			for key in opt_problem._variables.keys():
				if (opt_problem._variables[key].lower != -inf):
					g[i] = x[j] - xl[j]
					i += 1
				#end
				if (opt_problem._variables[key].upper != inf):
					g[i] = xu[j] - x[j]
					i += 1
				#end
				j += 1
			#end
			
			return f,g
		
		
		
		# Variables Handling
		nvar = len(opt_problem._variables.keys())
		xl = []
		xu = []
		xx = []
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'c'):
				xl.append(opt_problem._variables[key].lower)
				xu.append(opt_problem._variables[key].upper)
				xx.append(opt_problem._variables[key].value)
			elif (opt_problem._variables[key].type == 'i'):
				raise IOError('COBYLA cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('COBYLA cannot handle discrete design variables')
			#end
		#end
		xl = numpy.array(xl)
		xu = numpy.array(xu)
		xx = numpy.array(xx)
		
		# Variables Groups Handling 
		if opt_problem.use_groups:
			group_ids = {}
			k = 0
			for key in opt_problem._vargroups.keys():
				group_len = len(opt_problem._vargroups[key]['ids'])
				group_ids[opt_problem._vargroups[key]['name']] = [k,k+group_len]
				k += group_len
			#end
		#end		
		
		# Constraints Handling
		ncon = len(opt_problem._constraints.keys())
		#neqc = 0
		gg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					raise IOError('COBYLA cannot handle equality constraints')
					#neqc += 1
				#end
				#gg.append(opt_problem._constraints[key].value)
				gg.append(opt_problem._constraints[key].upper)
			#end
		#end
		nadd = 0
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].lower != -inf): 
				gg.append(0)
				nadd += 1
			#end
			if (opt_problem._variables[key].upper != inf): 
				gg.append(0)
				nadd += 1
			#end
		#end
		gg = numpy.array(gg,numpy.float)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff,numpy.float)
		
		
		# Setup argument list values
		n = numpy.array([nvar], numpy.int)
		m = numpy.array([ncon+nadd], numpy.int)
		rhobeg = numpy.array([self.options['RHOBEG'][1]], numpy.float)
		rhoend = numpy.array([self.options['RHOEND'][1]], numpy.float)
		if (myrank == 0):
			if (self.options['IPRINT'][1]>=0):
				iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
			else:
				raise IOError('Incorrect Output Level Setting')
			#end
		else:
			iprint = numpy.array([0], numpy.int)
		#end
		maxfun = numpy.array([self.options['MAXFUN'][1]], numpy.int)
		w = numpy.zeros([n*(3*n+2*m+11)+4*m+6], numpy.float)
		iact = numpy.zeros([m+1], numpy.intc)
		ifail = numpy.array([0], numpy.int)
		nfvals = numpy.array([0], numpy.int)
		iout = numpy.array([self.options['IOUT'][1]], numpy.int)
		ifile = self.options['IFILE'][1]
		if (iprint > 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		
		
		# Run COBYLA
		t0 = time.time()
		cobyla.cobyla(calcfc,n,m,xx,rhobeg,rhoend,iprint,maxfun,w,iact,
			ifail,nfvals,iout,ifile,ff,gg)
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
			cobyla.closeunit(self.options['IOUT'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = ifail[0]
		sol_inform['text'] = self.getInform(ifail[0])
		
		if store_sol:
			
			sol_name = 'COBYLA Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = nfvals[0]
			
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
			
			if ncon > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = -gg[i]
					i += 1
					if (i >= ncon):
						break
					#end
				#end
			else:
				sol_cons = {}
			#end
			
			sol_lambda = {}
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, myrank=myrank,
				arguments=args, **kwargs)
			
		#end
		
		return ff, xx, sol_inform
		
		
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 17, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getOption(self, name):
		
		'''
		Get Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 17, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getInform(self, infocode):
		
		'''
		Get Optimizer Result Information (Optimizer Specific Routine)
		
		Keyword arguments:
		-----------------
		id -> STRING: Option Name
		
		Documentation last updated:  May. 17, 2008 - Ruben E. Perez
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
			cobyla.pyflush(self.options['IOUT'][1])	
		#end
	


#==============================================================================
# COBYLA Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test COBYLA
	print 'Testing ...'
	COBYLA = COBYLA()
	print COBYLA
	
