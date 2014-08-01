#!/usr/bin/env python
'''
pyPSQP - A Python pyOpt interface to PSQP. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 31/07/2014 21:00$


Tested on:
---------
Win32 with g77
Linux with pathf95
Linux with gfortran
Linux with g77

Developers:
-----------
- Dr. Ruben E. Perez (RP)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2010)
	v. 1.1	- Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# PSQP Library
# =============================================================================
try:
	import psqp
except:
	raise ImportError('PSQP shared library failed to import')
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
# PSQP Optimizer Class
# =============================================================================
class PSQP(Optimizer):
	
	'''
	PSQP Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		PSQP Optimizer Class Initialization
		
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
		name = 'PSQP'
		category = 'Local Optimizer'
		def_opts = {
		'XMAX':[float,1e16],  		# Maximum Stepsize
		'TOLX':[float,1e-16],  		# Variable Change Tolerance
		'TOLC':[float,1e-6],  		# Constraint Violation Tolerance
		'TOLG':[float,1e-6],  		# Lagrangian Gradient Tolerance
		'RPF':[float,1e-4],  		# Penalty Coefficient
		'MIT':[int,1000],  			# Maximum Number of Iterations
		'MFV':[int,2000],  			# Maximum Number of Function Evaluations
		'MET':[int,2],  			# Variable Metric Update (1 - BFGS, 2 - Hoshino)
		'MEC':[int,2],  			# Negative Curvature Correction (1 - None, 2 - Powell's Correction)	
		'IPRINT':[int,2],			# Output Level (0 - None, 1 - Final, 2 - Iter)
		'IOUT':[int,6],     		# Output Unit Number
		'IFILE':[str,'PSQP.out'],	# Output File Name
		}
		informs = {
		1 : 'Change in design variable was less than or equal to tolerance',
		2 : 'Change in objective function was less than or equal to tolerance',
		3 : 'Objective function less than or equal to tolerance',
		4 : 'Maximum constraint value is less than or equal to tolerance',
		11 : 'Maximum number of iterations exceeded',
		12 : 'Maximum number of function evaluations exceeded',
		13 : 'Maximum number of gradient evaluations exceeded',
		-6 : 'Termination criterion not satisfied, but obtained point is acceptable',
		#<0 : 'Method failed',
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
		
		Documentation last updated:  February. 2, 2011 - Ruben E. Perez
		'''
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pyPSQP - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyPSQP: Parallel objective Function Analysis requires mpi4py'
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
		# PSQP - Objective/Constraint Values Storage 
		#======================================================================
		def eval(x):
			
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
			f = []
			g = []
			if (myrank == 0):
				if self.h_start:
					[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						[f,g,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
					#end
				#end
			#end
			
			if self.pll:
				self.h_start = Bcast(self.h_start,root=0)
			#end
			if self.h_start and self.pll:
				[f,g,fail] = Bcast([f,g,fail],root=0)
			elif not self.h_start:	
				[f,g,fail] = opt_problem.obj_fun(xn, *args, **kwargs)
			#end
			
			# Store History
			if (myrank == 0):
				if self.sto_hst:
					log_file.write(x,'x')
					log_file.write(f,'obj')
					log_file.write(g,'con')
					log_file.write(fail,'fail')
				#end
			#end
			
			# Objective Assigment
			if isinstance(f,float):
				f = [f]
			#end
			for i in xrange(len(opt_problem._objectives.keys())):
				if isinstance(f[i],complex):
					ff[i] = f[i].astype(float)
				else:
					ff[i] = f[i]
				#end
			#end
			
			# Constraints Assigment
			i = 0
			for j in xrange(len(opt_problem._constraints.keys())):
				if isinstance(g[j],complex):
					gg[i] = g[j].astype(float)
				else:
					gg[i] = g[j]
				#end
				i += 1
			#end
			
			# Gradients
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
				dff,dgg = gradient.getGrad(x, group_ids, f, g, *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')
			#end
			
			# Store
			self.stored_data['x'] = copy.copy(x)
			self.stored_data['f'] = copy.copy(ff)
			self.stored_data['g'] = copy.copy(gg)
			self.stored_data['df'] = copy.copy(dff)
			self.stored_data['dg'] = copy.copy(dgg)			
			
			return
		
		
		#======================================================================
		# PSQP - Objective Values Function
		#======================================================================
		def pobj(n,x,f):
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			ff = self.stored_data['f']
			
			return ff[0]
		
		
		#======================================================================
		# PSQP - Constraint Values Function
		#======================================================================
		def pcon(n,k,x,g):
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			gg = self.stored_data['g']
			
			return gg[k-1]
		
		
		#======================================================================
		# PSQP - Objective Gradients Function
		#======================================================================
		def pdobj(n,x,df):
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			df = self.stored_data['df']
			
			return df[0]
		
		
		#======================================================================
		# PSQP - Constraint Gradients Function
		#======================================================================
		def pdcon(n,k,x,dg):
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			dg = self.stored_data['dg']
			
			return dg[k-1]
		
		
		
		# Variables Handling
		nvar = len(opt_problem._variables.keys())
		xl = []
		xu = []
		xi = []
		xx = []
		for key in opt_problem._variables.keys():
			xl.append(opt_problem._variables[key].lower)
			xu.append(opt_problem._variables[key].upper)
			xi.append(3)
			xx.append(opt_problem._variables[key].value)
		#end
		xl = numpy.array(xl)
		xu = numpy.array(xu)
		xi = numpy.array(xi)
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
		if ncon > 0:
			gi = []
			gg = []
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					gi.append(5)
				elif opt_problem._constraints[key].type == 'i':
					gi.append(2)
				#end
				gg.append(opt_problem._constraints[key].value)
			#end
			gg.append(0.0)
			gl = numpy.zeros([ncon], numpy.float)
			gu = numpy.zeros([ncon], numpy.float)
			gi = numpy.array(gi, numpy.float)
			gg = numpy.array(gg, numpy.float)
		else:
			gl = numpy.array([0], numpy.float)
			gu = numpy.array([0], numpy.float)
			gi = numpy.array([0], numpy.float)
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
		nf = numpy.array([nvar], numpy.int)
		nc = numpy.array([ncon], numpy.int)
		mit = numpy.array([self.options['MIT'][1]], numpy.int)
		mfv = numpy.array([self.options['MFV'][1]], numpy.int)
		met = numpy.array([self.options['MET'][1]], numpy.int)
		mec = numpy.array([self.options['MEC'][1]], numpy.int)
		xmax = numpy.array([self.options['XMAX'][1]], numpy.float)
		tolx = numpy.array([self.options['TOLX'][1]], numpy.float)
		tolc = numpy.array([self.options['TOLC'][1]], numpy.float)
		tolg = numpy.array([self.options['TOLG'][1]], numpy.float)
		rpf = numpy.array([self.options['RPF'][1]], numpy.float)
		gmax = numpy.array([0], numpy.float)		
		cmax = numpy.array([0], numpy.float)
		if (myrank == 0):
			if (self.options['IPRINT'][1] <= 2):
				iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
			else:
				raise IOError('Incorrect Output Level Setting')
			#end
		else:
			iprint = numpy.array([0], numpy.int)
		#end
		iout = numpy.array([self.options['IOUT'][1]], numpy.int)
		ifile = self.options['IFILE'][1]
		if (myrank == 0):
			if (iprint != 0):
				if os.path.isfile(ifile):
					os.remove(ifile)
				#end
			#end
		#end
		iterm = numpy.array([0], numpy.int)
		
		
		# Storage Arrays 
		self.stored_data = {}
		self.stored_data['x'] = {}  #numpy.zeros([nvar],float)
		self.stored_data['f'] = {}  #numpy.zeros([nobj],float)
		self.stored_data['g'] = {}  #numpy.zeros([ncon],float)
		self.stored_data['df'] = {} #numpy.zeros([nvar],float)
		self.stored_data['dg'] = {} #numpy.zeros([ncon,nvar],float) 
		
		
		# Run PSQP
		t0 = time.time()
		psqp.psqp_wrap(nf,nc,xx,xi,xl,xu,gg,gi,gl,gu,mit,mfv,
			met,mec,xmax,tolx,tolc,tolg,rpf,ff,gmax,cmax,
			iprint,iout,ifile,iterm,pobj,pdobj,pcon,pdcon)
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
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = iterm[0]
		sol_inform['text'] = self.getInform(iterm[0])
		
		if store_sol:
			
			sol_name = 'PSQP Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = psqp.stat.nfv + psqp.stat.nfg*nvar
			
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
		iprint = self.options['IPRINT'][1]
		if (iprint > 0):
			psqp.pyflush(self.options['IOUT'][1])	
		#end
	


#==============================================================================
# PSQP Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test PSQP
	print 'Testing ...'
	psqp = PSQP()
	print psqp
	
