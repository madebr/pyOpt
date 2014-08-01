#!/usr/bin/env python
'''
pyMMA - A Python pyOpt interface to MMA. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.4   $Date: 21/06/2010 21:00$


Tested on:
---------
Linux with g77
Linux with gfortran
Linux with pathf95
Win32 with g77
Mac with g95

Developers:
-----------
- Mr. Andrew Lambe (AL)
- Dr. Ruben E. Perez (RP)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0	- Initial Class Creation (AL, 2005)
	v. 1.1	- Extensive Functionality Updates (RP, 2008)
	v. 1.2	- Migrate to pyOpt Framework (RP, 2008)
	v. 1.3	- History support (PJ,RP, 2010)
	v. 1.4	- Gradient Class Support (PJ,RP, 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- add unconstrained problems support
'''

# =============================================================================
# MMA Library
# =============================================================================
try:
	import mma
except:
	raise ImportError('MMA shared library failed to import')
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
# MMA Optimizer Class
# =============================================================================
class MMA(Optimizer):
	
	'''
	MMA Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		MMA Optimizer Class Initialization
		
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
		name = 'MMA'
		category = 'Local Optimizer'
		def_opts = {
		# MMA Options
		'MAXIT':[int,1000],     	# Maximum Iterations
		'GEPS':[float,1e-6],    	# Dual Objective Gradient Tolerance
		'DABOBJ':[float,1e-6],  	# 
		'DELOBJ':[float,1e-6],  	# 
		'ITRM':[int,2],         	# 
		'IPRINT':[int,1],       	# Output Level (<0 - None, 0 - Screen, 1 - File)
		'IOUT':[int,6],         	# Output Unit Number
		'IFILE':[str,'MMA.out'],	# Output File Name
		}
		informs = {
		0 : 'The optimality conditions are satisfied.', 
		1 : 'The algorithm has been stopped after MAXIT iterations.',
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
			raise NotImplementedError("pyMMA - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyMMA: Parallel objective Function Analysis requires mpi4py'
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
		# MMA - Objective/Constraint Values and Gradients Function
		#======================================================================
		def func(m,n,xval,f0val,df0dx,fval,dfdx):
			
			# Variables Groups Handling
			if opt_problem.use_groups:
				xg = {}
				for group in group_ids.keys():
					if (group_ids[group][1]-group_ids[group][0] == 1):
						xg[group] = xval[group_ids[group][0]]
					else:
						xg[group] = xval[group_ids[group][0]:group_ids[group][1]]
					#end
				#end
				xn = xg
			else:
				xn = xval
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
					log_file.write(xval,'x')
					log_file.write(f,'obj')
					log_file.write(g,'con')
					log_file.write(fail,'fail')
				#end
			#end
			
			# Gradients
			if self.h_start:
				df = []
				dg = []
				if (myrank == 0):
					[vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						df = vals['grad_obj'][0].reshape((len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())))
						dg = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))	
					#end
				#end
				if self.pll:
					self.h_start = Bcast(self.h_start,root=0)
				#end
				if self.h_start and self.pll:
					[df,dg] = Bcast([df,dg],root=0)
				#end
			#end
			
			if not self.h_start:
				
				#
				df,dg = gradient.getGrad(xval, group_ids, [f], g, *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(df,'grad_obj')
				log_file.write(dg,'grad_con')
			#end
			
			# Objective Assigment
			if isinstance(f,complex):
				f0val = f.astype(float)
			else:
				f0val = f
			#end
			
			# Constraints Assigment
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(g[i],complex):
					fval[i] = g[i].astype(float)
				else:
					fval[i] = g[i]
				#end
			#end
			
			# Gradients Assigment
			k = 0
			for i in xrange(len(opt_problem._variables.keys())):
				if isinstance(df[0,i],complex):
					df0dx[i] = df[0,i].astype(float)
				else:
					df0dx[i] = df[0,i]
				#end
				for jj in xrange(len(opt_problem._constraints.keys())):
					if isinstance(dg[jj,i],complex):
						dfdx[k] = dg[jj,i].astype(float)
					else:
						dfdx[k] = dg[jj,i]
					#end
					k += 1
				#end
			#end
			
			return f0val,df0dx,fval,dfdx
		
		
		
		# Variables Handling
		n = len(opt_problem._variables.keys())
		xmin = []
		xmax = []
		xval = []
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'c'):
				xmin.append(opt_problem._variables[key].lower)
				xmax.append(opt_problem._variables[key].upper)
				xval.append(opt_problem._variables[key].value)
			elif (opt_problem._variables[key].type == 'i'):
				raise IOError('MMA cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('MMA cannot handle discrete design variables')
			#end
		#end
		xmin = numpy.array(xmin)
		xmax = numpy.array(xmax)
		xval = numpy.array(xval)
		
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
		neqc = 0
		#fval = []
		fmax = []
		if m > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					raise IOError('MMA cannot handle equality constraints')
					#neqc += 1
				#end
				#fval.append(opt_problem._constraints[key].value)
				fmax.append(opt_problem._constraints[key].upper)
			#end
		else:
			raise IOError('MMA support for unconstrained problems not implemented yet')
		#end
		#fval = numpy.array(fval)
		fmax = numpy.array(fmax)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		f0val = []
		for key in opt_problem._objectives.keys():
			f0val.append(opt_problem._objectives[key].value)
		#end
		f0val = numpy.array(f0val)
		
		
		# Setup argument list values
		
		xmma = numpy.zeros([n], numpy.float)
		
		# Space used internally by the program 
		# for the asymptotes (xlow and xupp) and 
		# computed bounds on x (alpha and beta)
		xlow = numpy.zeros([n], numpy.float)
		xupp = numpy.zeros([n], numpy.float)
		alfa = numpy.zeros([n], numpy.float)
		beta = numpy.zeros([n], numpy.float)
		
		# The objective and constraint function 
		# values and space for the gradients
		fval = numpy.zeros([m], numpy.float)
		df0dx = numpy.zeros([n], numpy.float)
		dfdx = numpy.zeros([m*n], numpy.float)
		
		# Space for the coefficients and artificial 
		# variables to be computed (set to default values)
		p = numpy.zeros([m*n], numpy.float)
		q = numpy.zeros([m*n], numpy.float)
		p0 = numpy.zeros([n], numpy.float)
		q0 = numpy.zeros([n], numpy.float)
		b = numpy.zeros([m], numpy.float)
		y = numpy.zeros([m], numpy.float)
		z = numpy.array([0.], numpy.float)
		a = numpy.zeros([m], numpy.float)
		c = 10000*numpy.ones([m], numpy.float)
		
		# Space for the Lagrange multipliers (ulam) 
		# the gradient of the dual objective function,
		# search direction, and Hessian of the dual objective
		ulam = numpy.ones([m], numpy.float)
		gradf = numpy.zeros([m], numpy.float)
		dsrch = numpy.zeros([m], numpy.float)
		hessf = numpy.zeros([m*(m+1)/2], numpy.float)
		
		# Specify that all variables are free to move
		iyfree = numpy.ones([m], numpy.int)
		
		# 
		iter = numpy.array([0], numpy.int)
		maxit = numpy.array([self.options['MAXIT'][1]], numpy.int)
		geps = numpy.array([self.options['GEPS'][1]], numpy.float)
		dabobj = numpy.array([self.options['DABOBJ'][1]], numpy.float)
		delobj = numpy.array([self.options['DELOBJ'][1]], numpy.float)
		itrm = numpy.array([self.options['ITRM'][1]], numpy.int)
		inform = numpy.array([0], numpy.int)
		if (myrank == 0):
			iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
		else:
			iprint = numpy.array([0], numpy.int)
		iout = numpy.array([self.options['IOUT'][1]], numpy.int)
		ifile = self.options['IFILE'][1]
		if (myrank == 0):
			if (iprint >= 0):
				if os.path.isfile(ifile):
					os.remove(ifile)
				#end
			#end
		#end
		
		#
		nfunc = numpy.array([0], numpy.int)
		
		# Run MMA
		t0 = time.time()
		mma.mma(n,m,iter,maxit,geps,dabobj,delobj,itrm,inform,
			xval,xmma,xmin,xmax,xlow,xupp,alfa,beta,f0val,fval,
			fmax,df0dx,dfdx,p,q,p0,q0,b,y,z,a,c,ulam,gradf,
			dsrch,hessf,iyfree,iprint,iout,ifile,nfunc,func)
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
			mma.closeunit(self.options['IOUT'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = inform[0]
		sol_inform['text'] = self.getInform(inform[0])
		
		if store_sol:
			
			sol_name = 'MMA Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = nfunc[0]
			
			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = xmma[i]
				i += 1
			#end
			
			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = f0val[i]
				i += 1
			#end
			
			if m > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = fval[i]
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
		
		return f0val, xmma, sol_inform
		
		
		
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
			mma.pyflush(self.options['IOUT'][1])
		#end
	


#==============================================================================
# MMA Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test MMA
	print 'Testing ...'
	mma = MMA()
	print mma
	
