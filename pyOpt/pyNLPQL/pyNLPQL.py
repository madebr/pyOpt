#!/usr/bin/env python
'''
pyNLPQL - A Python pyOpt interface to NLPQL. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.5   $Date: 21/06/2010 21:00$


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
- Mr. C.A.(Sandy) Mader (SM)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0 	- Initial Class Creation (RP, 2007)
	v. 1.1 	- Wrapper (f2py) Fixes (SM, 2007)
			- Interface Fixes (RP, 2007)
	v. 1.2 	- Migrate to pyOpt Framework (RP, 2008)
			- Bug fix for CS (PJ, 2008)
	v. 1.3  - User-Provided Sensitivities Support (PJ,RP, 2008)
	v. 1.4  - History support (PJ,RP, 2010)
	v. 1.5	- Gradient Class Support (PJ,RP, 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# NLPQL Library
# =============================================================================
try:
	import nlpql
except:
	raise ImportError('NLPQL shared library failed to import')
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
# NLPQL Optimizer Class
# =============================================================================
class NLPQL(Optimizer):
	
	'''
	NLPQL Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		NLPQL Optimizer Class Initialization
		
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
		name = 'NLPQL'
		category = 'Local Optimizer'
		def_opts = {
		# NLPQL Options
		'Accurancy':[float,1e-6],   # Convergence Accurancy
		'ScaleBound':[float,1e30],  # 
		'maxFun':[int,20],          # Maximum Number of Function Calls During Line Search
		'maxIt':[int,500],          # Maximum Number of Iterations
		'iPrint':[int,2],           # Output Level (0 - None, 1 - Final, 2 - Major, 3 - Major/Minor, 4 - Full)
		'mode':[int,0],             # NLPQL Mode (0 - Normal Execution, 1 to 18 - See Manual)
		'iout':[int,6],             # Output Unit Number
		'lmerit':[bool,True],       # Merit Function Type (True - L2 Augmented Penalty, False - L1 Penalty)
		'lql':[bool,False],         # QP Subproblem Solver (True - Quasi-Newton, False - Cholesky)
		'iFile':[str,'NLPQL.out'],	# Output File Name
		}
		informs = {
		-2 : 'Compute gradient values w.r.t. the variables stored in' \
			' first column of X, and store them in DF and DG.' \
			' Only derivatives for active constraints ACTIVE(J)=.TRUE. need to be computed.',
		-1 : 'Compute objective fn and all constraint values subject' \
			'the variables found in the first L columns of X, and store them in F and G.',
		0 : 'The optimality conditions are satisfied.', 
		1 : ' The algorithm has been stopped after MAXIT iterations.',
		2 : ' The algorithm computed an uphill search direction.',
		3 : ' Underflow occurred when determining a new approximation matrix' \
			'for the Hessian of the Lagrangian.',
		4 : 'The line search could not be terminated successfully.', 
		5 : 'Length of a working array is too short.' \
			' More detailed error information is obtained with IPRINT>0',
		6 : 'There are false dimensions, for example M>MMAX, N>=NMAX, or MNN2<>M+N+N+2.',
		7 : 'The search direction is close to zero, but the current iterate is still infeasible.',
		8 : 'The starting point violates a lower or upper bound.',
		9 : 'Wrong input parameter, i.e., MODE, LDL decomposition in D and C' \
			' (in case of MODE=1), IPRINT, IOUT',
		10 : 'Internal inconsistency of the quadratic subproblem, division by zero.',
		100 : 'The solution of the quadratic programming subproblem has been' \
			' terminated with an error message and IFAIL is set to IFQL+100,' \
			' where IFQL denotes the index of an inconsistent constraint.',
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
			raise NotImplementedError("pyNLPQL - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyNLPQL: Parallel objective Function Analysis requires mpi4py'
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
		def_fname = self.options['iFile'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		#
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)
		
		
		#======================================================================
		# NLPQL - Objective/Constraint Values Function (Real Valued) 
		#======================================================================
		def nlfunc(m,me,mmax,n,f,g,x,active):
			
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
			
			# Constraints Assigment (negative gg as nlpql uses g(x) >= 0)
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(gg[i],complex):
					g[i] = -gg[i].astype(float)
				else:
					g[i] = -gg[i]
				#end
			#end
			
			return f,g
		
		
		#======================================================================
		# NLPQL - Objective/Constraint Gradients Function
		#======================================================================
		def nlgrad(m,me,mmax,n,f,g,df,dg,x,active,wa):
			
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
				dff,dgg = gradient.getGrad(x, group_ids, [f], -g[0:len(opt_problem._constraints.keys())], *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')
			#end
			
			# Gradient Assignment
			for i in xrange(len(opt_problem._variables.keys())):
				df[i] = dff[0,i]
				for j in xrange(len(opt_problem._constraints.keys())):
					dg[j,i] = -dgg[j,i]
				#end
			#end
			
			return df,dg
		
		
		
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
				raise IOError('NLPQL cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('NLPQL cannot handle discrete design variables')
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
		ncon = len(opt_problem._constraints.keys())
		neqc = 0
		#gg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					neqc += 1
				#end
				#gg.append(opt_problem._constraints[key].value)
			#end
		#end
		#gg = numpy.array(gg, numpy.float)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff, numpy.float)
		
		
		# Setup argument list values
		mm = numpy.array([ncon], numpy.int)
		me = numpy.array([neqc], numpy.int)
		mmx = 200
		if (ncon >= mmx):
			mmx = ncon + 1
		#end
		mmax = numpy.array([mmx], numpy.int)
		nn = numpy.array([nvar], numpy.int)
		nmx = 200
		if (nvar >= nmx):
			nmx = nvar + 1
		#end
		nmax = numpy.array([nmx], numpy.int)
		mnn2 = numpy.array([mm[0]+nn[0]+nn[0]+2], numpy.int)
		#xx = _ConcatenateVector(self.variables, 'value')
		#ff = self.objective.value
		gg = numpy.zeros([mmax], numpy.float)
		df = numpy.zeros([nmax], numpy.float)
		dg = numpy.zeros([mmax,nmax], numpy.float)
		uu = numpy.zeros([mnn2], numpy.float)
		#xl = _ConcatenateVector(self.variables, 'lower')
		#xu = _ConcatenateVector(self.variables, 'upper')
		cc = numpy.zeros([nmax,nmax], numpy.float)
		dd = numpy.zeros([nmax], numpy.float)
		acc = numpy.array([self.options['Accurancy'][1]], numpy.float)
		scbou = numpy.array([self.options['ScaleBound'][1]], numpy.float)
		maxfun = numpy.array([self.options['maxFun'][1]], numpy.int)
		maxit = numpy.array([self.options['maxIt'][1]], numpy.int)
		if (myrank == 0):
			if (self.options['iPrint'][1]>=0 and self.options['iPrint'][1]<=4):
				iprint = numpy.array([self.options['iPrint'][1]], numpy.int)
			else:
				raise IOError('Incorrect Output Level Setting')
			#end
		else:
			iprint = numpy.array([0], numpy.int)
		#end
		if (self.options['mode'][1]>=0 and self.options['mode'][1]<=18):
			mode = self.options['mode'][1]
		else:
			raise IOError('Incorrect Mode Setting')
		#end
		iout = self.options['iout'][1]
		ifile = self.options['iFile'][1]
		if (iprint > 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		ifail = numpy.array([0], numpy.int)
		lwa0 = 100000
		lwa1 = 4*mmax + 4*ncon + 19*nvar + 55
		lwa2 = mmax*nvar + 4*mmax + 4*ncon + 18*nvar + 55
		lwa3 = 3/2*(nvar + 1)*(nvar + 1) + 10*nvar + 2*ncon + 10
		lwaM = max([lwa0,lwa1,lwa2,lwa3]) + lwa3
		lwa = numpy.array([lwaM], numpy.int)
		wa = numpy.zeros([lwa], numpy.float)
		lkwa = numpy.array([mmx+2*nmx+20], numpy.int)
		kwa = numpy.zeros([lkwa], numpy.intc)
		lactiv = numpy.array([2*mmx+15], numpy.int)
		active = numpy.zeros([lactiv], numpy.bool)
		lmerit = numpy.array([self.options['lmerit'][1]], numpy.bool)
		lql = numpy.array([self.options['lql'][1]], numpy.bool)
		fmp = numpy.array([eps], numpy.float)
		
		
		# Run NLPQL
		t0 = time.time()
		nlpql.nlpql1(mm,me,mmax,nn,nmax,mnn2,xx,ff,gg,df,dg,uu,xl,xu,cc,
			dd,acc,scbou,maxfun,maxit,iprint,mode,iout,ifile,ifail,
			wa,lwa,kwa,lkwa,active,lactiv,lmerit,lql,fmp,nlfunc,nlgrad)
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
			nlpql.closeunit(self.options['iout'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = ifail[0]
		sol_inform['text'] = self.getInform(ifail[0])
		
		if store_sol:
			
			sol_name = 'NLPQL Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = kwa[0] + kwa[1]*nvar
			
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
				#end
			else:
				sol_cons = {}
			#end
			
			if ncon > 0:
				sol_lambda = numpy.zeros(ncon,float)
				for i in xrange(ncon):
					sol_lambda[i] = uu[i]
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
		
		# 
		if (infocode <= 10):
			return self.informs[infocode]
		else:
			return self.informs[100]
		#end
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		'''
		
		# 
		iPrint = self.options['iPrint'][1]
		if (iPrint > 0):
			nlpql.pyflush(self.options['iout'][1])
		#end
	


#==============================================================================
# NLPQL Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test NLPQL
	print 'Testing ...'
	nlpql = NLPQL()
	print nlpql
	
