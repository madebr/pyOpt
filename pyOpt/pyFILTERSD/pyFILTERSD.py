#!/usr/bin/env python
'''
pyFILTERSD - A Python pyOpt interface to filterSD. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with gfortran
Win32 with gfortran

Developers:
-----------
- Dr. Ruben E. Perez (RP)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2012)
	v. 1.1	- Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# FILTERSD Library
# =============================================================================
try:
	import filtersd
except:
	raise ImportError('FILTERSD shared library failed to import')
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
# FILTERSD Optimizer Class
# =============================================================================
class FILTERSD(Optimizer):
	
	'''
	FILTERSD Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		FILTERSD Optimizer Class Initialization
		
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
		name = 'FILTERSD'
		category = 'Local Optimizer'
		def_opts = {
		'rho':[float,100.0],			# initial trust region radius
		'htol':[float,1e-6],			# tolerance allowed in sum h of constraint feasibilities
		'rgtol':[float,1e-5],			# tolerance allowed in reduced gradient l2 norm
		'maxit':[int,1000],				# maximum number of major iterations allowed
		'maxgr':[int,1e5],				# upper limit on the number of gradient calls
		'ubd':[float,1e5],				# upper bound on the allowed constraint violation
		'dchk':[int,0],					# derivative check flag (0 - no check, 1 - check)
		'dtol':[float,1e-8],			# derivative check tolerance
		'iprint':[int,1],				# verbosity of printing (0 - none, 1 - Iter, 2 - Debug)
		'iout':[int,6],     			# Output Unit Number
		'ifile':[str,'FILTERSD.out'],	# Output File Name
		}
		informs = {
		-1 : 'ws not large enough',
		-2 : 'lws not large enough',
		-3 : 'inconsistency during derivative check',
		0 : 'successful run',
		1 : 'unbounded NLP (f <= fmin at an htol-feasible point)',
		2 : 'bounds on x are inconsistent',
		3 : 'local minimum of feasibility problem and h > htol, (nonlinear constraints are locally inconsistent)',
		4 : 'initial point x has h > ubd (reset ubd or x and re-enter)',
		5 : 'maxit major iterations have been carried out',
		6 : 'termination with rho <= htol',
		7 : 'not enough workspace in ws or lws (see message)',
		8 : 'insufficient space for filter (increase mxf and re-enter)',
		9 : 'unexpected fail in LCP solver',
		10 : 'unexpected fail in LCP solver',
		}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		
	def __solve__(self, opt_problem={}, sens_type='FD', store_sol=True, store_hst=False, hot_start=False, disp_opts=False, sens_mode='', sens_step={}, *args, **kwargs):
		
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
		
		Documentation last updated:  February. 2, 2013 - Ruben E. Perez
		'''
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pyFILTERSD - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyFILTERSD: Parallel objective Function Analysis requires mpi4py'
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
		def_fname = self.options['ifile'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		# 
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)
		
		
		#======================================================================
		# filterSD - Objective/Constraint Values Function (Real Valued) 
		#======================================================================
		def functions(n,m,x,f,g,user,iuser):
			
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
			
			# Store
			self.stored_data['x'] = copy.copy(x)
			self.stored_data['f'] = copy.copy(ff)
			self.stored_data['g'] = copy.copy(gg)
			
			# Objective Assigment
			if isinstance(ff,complex):
				f = ff.astype(float)
			else:
				f = ff
			#end
			
			# Constraints Assigment
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(gg[i],complex):
					g[i] = gg[i].astype(float)
				else:
					g[i] = gg[i]
				#end
			#end
			
			return f,g
		
		
		#======================================================================
		# filterSD - Objective/Constraint Gradients Function
		#======================================================================
		def gradients(n,m,x,a,user,iuser):
			
			if ((self.stored_data['x'] != x).any()):
				f,g = functions(n,m,x,[],[],[],[])
			else:
				f = self.stored_data['f']
				g = self.stored_data['g']
			#end
			
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
				dff,dgg = gradient.getGrad(x, group_ids, [f], g[0:len(opt_problem._constraints.keys())], *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')
			#end
			
			# Gradient Assignment
			for i in xrange(len(opt_problem._variables.keys())):
				a[i,0] = dff[0,i]
			#end
			for i in xrange(len(opt_problem._variables.keys())):
				for j in xrange(len(opt_problem._constraints.keys())):
					a[i,j+1] = dgg[j,i]
				#end
			#end
			
			return a
		
		
		# Variables Handling
		nvar = len(opt_problem._variables.keys())
		xl = []
		xu = []
		xx = []
		for key in opt_problem._variables.keys():
			xl.append(opt_problem._variables[key].lower)
			xu.append(opt_problem._variables[key].upper)
			xx.append(opt_problem._variables[key].value)
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
		ncon = len(opt_problem._constraints.keys())
		gg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					raise IOError('FILTERSD cannot handle equality constraints')
				#end
				gg.append(opt_problem._constraints[key].value)
			#end
			gg = numpy.array(gg, numpy.float)
		else:
			ncon = 1
			gg = numpy.array([0], numpy.float)
		#end
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff, numpy.float)
		
		
		# Setup argument list values
		nn = numpy.array([nvar], numpy.int)
		mm = numpy.array([ncon], numpy.int)
		al = numpy.zeros([mm], numpy.float)
		ubd = numpy.array([self.options['ubd'][1]], numpy.float)
		rho = numpy.array([self.options['rho'][1]], numpy.float)
		htol = numpy.array([self.options['htol'][1]], numpy.float)
		rgtol = numpy.array([self.options['rgtol'][1]], numpy.float)
		maxit = numpy.array([self.options['maxit'][1]], numpy.int)
		maxgr = numpy.array([self.options['maxgr'][1]], numpy.int)
		dchk = numpy.array([self.options['dchk'][1]], numpy.int)
		dtol = numpy.array([self.options['dtol'][1]], numpy.float)
		if (myrank == 0):
			if (self.options['iprint'][1]>=0):
				iprint = numpy.array([self.options['iprint'][1]], numpy.int)
			else:
				raise IOError('Incorrect Output Level Setting')
			#end
		else:
			iprint = numpy.array([0], numpy.int)
		#end
		iout = numpy.array([self.options['iout'][1]], numpy.int)
		ifile = self.options['ifile'][1]
		if (iprint > 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		ifail = numpy.array([0], numpy.int)
		nfevs = numpy.array([0], numpy.int)
		ngevs = numpy.array([0], numpy.int)
		
		
		# Storage Arrays 
		self.stored_data = {}
		self.stored_data['x'] = {}
		self.stored_data['f'] = {}
		self.stored_data['g'] = {}
		
		
		# Run filterSD
		t0 = time.time()
		filtersd.filtersd_wrap(nn,mm,xx,xl,xu,al,ff,gg,inf,ubd,rho,htol,rgtol,maxit,maxgr,dchk,dtol,iprint,iout,ifile,ifail,nfevs,ngevs,functions,gradients)
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
			filtersd.closeunit(self.options['iout'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = ifail[0]
		sol_inform['text'] = self.getInform(ifail[0])
		
		if store_sol:
			
			sol_name = 'FILTERSD Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = nfevs[0] + ngevs[0]*nvar
			
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
					sol_cons[key].value = gg[i]
					i += 1
				#end
			else:
				sol_cons = {}
			#end
			
			if ncon > 0:
				sol_lambda = numpy.zeros(ncon,float)
				for i in xrange(ncon):
					sol_lambda[i] = al[i]
				#end
			else:
				sol_lambda = {}
			#end
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, Sensitivities=sens_type, 
				myrank=myrank, arguments=args, **kwargs)
			
			time.sleep(0)
		#end
		
		return ff, xx, sol_inform
		
		
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  November. 30, 2010 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getOption(self, name):
		
		'''
		Get Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  November. 30, 2010 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getInform(self, infocode):
		
		'''
		Get Optimizer Result Information (Optimizer Specific Routine)
		
		Keyword arguments:
		-----------------
		id -> STRING: Option Name
		
		Documentation last updated:  November. 30, 2010 - Ruben E. Perez
		'''
		
		return self.informs[infocode]
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  November. 30, 2010 - Ruben E. Perez
		'''
		
		# 
		iprint = self.options['iprint'][1]
		if (iprint > 0):
			filtersd.pyflush(self.options['iout'][1])	
		#end
	


#==============================================================================
# FILTERSD Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test FILTERSD
	print 'Testing ...'
	filtersd = FILTERSD()
	print filtersd
	
