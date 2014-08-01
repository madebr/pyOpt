#!/usr/bin/env python
'''
pyFSQP - A Python pyOpt interface to FSQP. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.5   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with g77
Win32 with g77
Linux with pathf95

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Andrew Lambe (AL)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2009)
	v. 1.1	- Integrate to pyOpt Framework (RP, 2009)
	v. 1.2  - Wrapper Callback Storage Support (AL,RP, 2009)
	v. 1.3	- History support (PJ,RP, 2010)
	v. 1.4  - Gradient Class Support (PJ,RP, 2010)
	v. 1.5  - Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- check possible ff,gg issue in eval
'''

# =============================================================================
# FSQP Library
# =============================================================================
try:
	import ffsqp
except:
	raise ImportError('FSQP shared library failed to import')
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
# FSQP Optimizer Class
# =============================================================================
class FSQP(Optimizer):
	
	'''
	FSQP Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		FSQP Optimizer Class Initialization
		
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
		name = 'FSQP'
		category = 'Local Optimizer'
		def_opts = {
		'mode':[int,100],			# FSQP Mode (See Manual)
		'iprint':[int,2],			# Output Level (0 - None, 1 - Final, 2 - Major, 3 - Major Details)
		'miter':[int,500],			# Maximum Number of Iterations
		'bigbnd':[float,1e10],		# Plus Infinity Value
		'epstol':[float,1e-8],		# Convergence Tolerance
		'epseqn':[float,0],			# Equality Constraints Tolerance
		'iout':[int,6],     		# Output Unit Number
		'ifile':[str,'FSQP.out'],	# Output File Name
		}
		informs = {
		0 : 'Normal termination of execution',
		1 : 'User-provided initial guess is infeasible for linear constraints, unable to generate a point satisfying all these constraints',
		2 : 'User-provided initial guess is infeasible for nonlinear inequality constraints and linear constraints, unable to generate a point satisfying all these constraints',
		3 : 'The maximum number of iterations has been reached before a solution is obtained',
		4 : 'The line search fails to find a new iterate',
		5 : 'Failure of the QP solver in attempting to construct d0, a more robust QP solver may succeed',
		6 : 'Failure of the QP solver in attempting to construct d1, a more robust QP solver may succeed',
		7 : 'Input data are not consistent, check print out error messages',
		8 : 'Two consecutive iterates are numerically equivalent before a stopping criterion is satisfied',
		9 : 'One of the penalty parameters exceeded bigbnd, the algorithm is having trouble satisfying a nonlinear equality constraint',
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
		nec = 0
		nic = 0
		for key in opt_problem._constraints.keys():
			if opt_problem._constraints[key].type == 'e':
				nec += 1
			if opt_problem._constraints[key].type == 'i':
				nic += 1
			#end
		#end
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pyFSQP - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyFSQP: Parallel objective Function Analysis requires mpi4py'
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
		# FSQP - Objective/Constraint Values Storage 
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
		# FSQP - Objective Values Function
		#======================================================================
		def obj(nparam,j,x,fj):
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			ff = self.stored_data['f']
			if (nobj == 1):
				fj = ff
			else:
				if isinstance(ff,list):
					ff = numpy.array(ff)
				#end
				fj = ff[j-1]
			#end
			
			return fj
		
		
		#======================================================================
		# FSQP - Constraint Values Function
		#======================================================================
		def cntr(nparam,j,x,gj):
			
			# for given j, assign to gj the value of the jth constraint evaluated at x
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end			
			
			gg = self.stored_data['g']
			if (j <= nic):
				jg = nec + (j-1)
			else:
				jg = (j-1) - nic
			#end
			gj = gg[jg]	
			
			return gj
		
		
		#======================================================================
		# FSQP - Objective Gradients Function
		#======================================================================
		def gradobj(nparam,j,x,gradfj,obj):
			
			# assign to gradfj the gradient of the jth objective function evaluated at x
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end			
			
			df = self.stored_data['df']
			for i in xrange(len(opt_problem._variables.keys())):
				gradfj[i] = df[j-1,i]
			#end
			
			return gradfj
		
		
		#======================================================================
		# FSQP - Constraint Gradients Function
		#======================================================================
		def gradcntr(nparam,j,x,gradgj,obj):
			
			# assign to gradgj the gradient of the jth constraint evaluated at x
			
			if ((self.stored_data['x'] != x).any()):
				eval(x)
			#end
			
			dg = self.stored_data['dg']
			if (j <= nic):
				jg = nec + (j-1)
			else:
				jg = (j-1) - nic
			#end
			gradgj = dg[jg]
			
			return gradgj
		
		
		
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
		neqc = 0
		gg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					neqc += 1
				#end
				gg.append(opt_problem._constraints[key].value)
			#end
			gg = numpy.array(gg, numpy.float)
		else:
			gg = numpy.array([0] ,numpy.float)
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
		nparam = numpy.array([nvar], numpy.int)
		nf = numpy.array([nobj], numpy.int)
		nineqn = numpy.array([ncon-neqc], numpy.int)
		nineq = nineqn
		neqn = numpy.array([neqc], numpy.int)
		neq = neqn
		mode = numpy.array([self.options['mode'][1]], numpy.int)
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
		miter = numpy.array([self.options['miter'][1]], numpy.int)
		inform = numpy.array([0], numpy.int)
		bigbnd = numpy.array([self.options['bigbnd'][1]], numpy.float)
		epstol = numpy.array([self.options['epstol'][1]], numpy.float)
		epsneq = numpy.array([self.options['epseqn'][1]], numpy.float)
		udelta = numpy.array([0], numpy.float)
		iwsizeM = 6*nvar + 8*max([1,ncon]) + 7*max([1,nobj]) + 30
		iwsize = numpy.array([iwsizeM], numpy.int)
		iw = numpy.zeros([iwsize], numpy.float)
		nwsizeM = 4*nvar**2 + 5*max([1,ncon])*nvar + 3*max([1,nobj])*nvar + 26*(nvar+max([1,nobj])) + 45*max([1,ncon]) + 100		
		nwsize = numpy.array([nwsizeM], numpy.int)
		w = numpy.zeros([nwsize], numpy.float)
		
		
		# Storage Arrays 
		self.stored_data = {}
		self.stored_data['x'] = {}  #numpy.zeros([nvar],float)
		self.stored_data['f'] = {}  #numpy.zeros([nobj],float)
		self.stored_data['g'] = {}  #numpy.zeros([ncon],float)
		self.stored_data['df'] = {} #numpy.zeros([nvar],float)
		self.stored_data['dg'] = {} #numpy.zeros([ncon,nvar],float) 
		
		
		# Run FSQP
		t0 = time.time()
		ffsqp.ffsqp(nparam,nf,nineqn,nineq,neqn,neq,mode,iprint,miter,
			inform,bigbnd,epstol,epsneq,udelta,xl,xu,xx,ff,gg,iw,iwsize,
			w,nwsize,obj,cntr,gradobj,gradcntr,iout,ifile)
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
			ffsqp.closeunit(self.options['iout'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = inform[0]
		sol_inform['text'] = self.getInform(inform[0])
		
		if store_sol:
			
			sol_name = 'FSQP Solution to ' + opt_problem.name
			
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
					sol_lambda[i] = w[nvar+i]
				#end
			else:
				sol_lambda = {}
			#end
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, Sensitivities=sens_type, 
				myrank=myrank, arguments=args, **kwargs)
			
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
		iprint = self.options['iprint'][1]
		if (iprint > 1):
			ffsqp.pyflush(self.options['iout'][1])	
		#end
	


#==============================================================================
# FSQP Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test FSQP
	print 'Testing ...'
	fsqp = FSQP()
	print fsqp
	
