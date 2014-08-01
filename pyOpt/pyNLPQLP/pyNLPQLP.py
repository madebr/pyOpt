#!/usr/bin/env python
'''
pyNLPQLP - A Python pyOpt interface to NLPQLP. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 18/12/2012 21:00$


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
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# NLPQLP Library
# =============================================================================
try:
	import nlpqlp
except:
	raise ImportError('NLPQLP shared library failed to import')
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
# NLPQLP Optimizer Class
# =============================================================================
class NLPQLP(Optimizer):
	
	'''
	NLPQLP Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		NLPQLP Optimizer Class Initialization
		
		**Keyword arguments:**
		
		- pll_type -> STR: Parallel Implementation (None, 'POA'-Parallel Objective Analysis), *Default* = None
		
		Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
		'''
		
		#
		if (pll_type == None):
			self.pll_type = None
		elif (pll_type.upper() == 'SPM'):
			self.pll_type = 'SPM'
		elif (pll_type.upper() == 'POA'):
			self.pll_type = 'POA'
		else:
			raise ValueError("pll_type must be either None,'SPM', 'DPM' or 'POA'")
		#end
		
		#
		name = 'NLPQLP'
		category = 'Local Optimizer'
		def_opts = {
		'ACC':[float,1e-8],         # Convergence Accurancy
		'ACCQP':[float,1e-12],      # QP Solver Convergence Accurancy
		'STPMIN':[float,1e-10],     # Minimum Step Length
		'MAXFUN':[int,10],          # Maximum Number of Function Calls During Line Search
		'MAXIT':[int,100],          # Maximum Number of Outer Iterations
		'RHOB':[float,0.0],         # BFGS-Update Matrix Initialization Parameter
		'IPRINT':[int,2],           # Output Level (0 - None, 1 - Final, 2 - Major, 3 - Major/Minor, 4 - Full)
		'MODE':[int,0],             # NLPQLP Mode
		'IOUT':[int,6],             # Output Unit Number
		'LQL':[bool,True],          # QP Subproblem Solver (True - Quasi-Newton, False - Cholesky)
		'IFILE':[str,'NLPQLP.out'], # Output File Name
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
		
		Documentation last updated:  February. 2, 2013 - Peter W. Jansen
		'''
		
		# 
		if ((self.pll_type != None) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pyNLPQLP - Current implementation only allows single level parallelization, either 'POA', 'SMP', 'DPM' or 'pgc'")
		#end
		
		if (self.pll_type != None) or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyNLPQLP: Parallel objective Function Analysis or gradient calculation requires mpi4py'
			#end
			comm = MPI.COMM_WORLD
			if (mpi4py.__version__[0] == '0'):
				Bcast = comm.Bcast
				Barrier = comm.Barrier
				Send = comm.Send
				Recv = comm.Recv
			elif (mpi4py.__version__[0] == '1'):
				Bcast = comm.bcast
				Barrier = comm.barrier
				Send = comm.send
				Recv = comm.recv
			#end
			self.pll = True
			if (self.pll_type == 'SPM'):
				nproc = comm.Get_size()
			else:
				nproc = 1
			#end
			self.myrank = comm.Get_rank()
		else:
			self.pll = False
			nproc = 1
			self.myrank = 0
		#end
		
		myrank = self.myrank
		
		# 
		def_fname = self.options['IFILE'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		# 
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)
		
		
		#======================================================================
		# NLPQLP - Objective/Constraint Values Function
		#======================================================================
		def nlfunc(l,nmax,mmax,x,lactive,active,f,g):
			
			#
			if (self.pll_type == 'SPM'):
				mxi = myrank
			else:
				mxi = 0
			#end
			
			# Variables Groups Handling
			if opt_problem.use_groups:
				xg = {}
				for group in group_ids.keys():
					if (group_ids[group][1]-group_ids[group][0] == 1):
						xg[group] = x[group_ids[group][0],mxi]
					else:
						xg[group] = x[group_ids[group][0]:group_ids[group][1],mxi]
					#end
				#end
				xn = xg
			else:
				xn = x[:-1,mxi]
			#end
			
			# Flush Output Files
			self.flushFiles()
				
			# Evaluate User Function
			fail = 0
			ff = []
			gg = []
			if (myrank == 0):
				if self.h_start:
					for proc in xrange(l):
						[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
						if hist_end:
							self.h_start = False
							hos_file.close()
						else:
							[ff,gg,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
							f[proc] = ff
							g[:,proc] = gg
						#end
					#end
				#end
			#end
				
			if self.pll:
				self.h_start = Bcast(self.h_start,root=0)
			#end
			if self.h_start and self.pll:
				[f,g] = Bcast([f,g],root=0)
			elif not self.h_start:	
				[ff,gg,fail] = opt_problem.obj_fun(xn, *args, **kwargs)
				
				# Objective Assigment
				if isinstance(ff,complex):
					f[mxi] = ff.astype(float)
				else:
					f[mxi] = ff
				#end
				
				# Constraints Assigment (negative gg as nlpqlp uses g(x) >= 0)
				for i in xrange(len(opt_problem._constraints.keys())):
					if isinstance(gg[i],complex):
						g[i,mxi] = -gg[i].astype(float)
					else:
						g[i,mxi] = -gg[i]
					#end
				#end
				
				if (self.pll_type == 'SPM'):
					send_buf = {}
					send_buf[myrank] = {'fi':f[mxi],'gi':g[:,mxi]}
					if myrank != 0:
						Send(send_buf,dest=0)
					else:
						p_results = []
						for proc in xrange(1,nproc):
							p_results.append(Recv(source=proc))
						#end
					#end
					
					if myrank == 0:
						for proc in xrange(nproc-1):
							for i in p_results[proc].keys():
								f[i] = p_results[proc][i]['fi']
								g[:,i] = p_results[proc][i]['gi']
							#end
						#end
					#end
					
					[f,g] = Bcast([f,g],root=0)
				#end
				
			#end
			
			# Store History
			if (myrank == 0):
				if self.sto_hst:
					for proc in xrange(l):
						log_file.write(x[:-1,proc],'x')
						log_file.write(f[proc],'obj')
						log_file.write(g[:,proc],'con')
						log_file.write(fail,'fail')
					#end
				#end
			#end

			
			return f,g
		
		
		#======================================================================
		# NLPQLP - Objective/Constraint Gradients Function
		#======================================================================
		def nlgrad(l,nmax,mmax,x,lactive,active,f,g,df,dg):
			
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
				dff,dgg = gradient.getGrad(x[:-1,0], group_ids, [f[0]], -g[0:len(opt_problem._constraints.keys()),0], *args, **kwargs)
				
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
		xl = numpy.zeros([max(2,nvar+1)], numpy.float)
		xu = numpy.zeros([max(2,nvar+1)], numpy.float)
		xx = numpy.zeros([max(2,nvar+1),nproc], numpy.float)
		i = 0
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'c'):
				xl[i] = opt_problem._variables[key].lower
				xu[i] = opt_problem._variables[key].upper
				for proc in xrange(nproc):
					xx[i,proc] = opt_problem._variables[key].value
				#end
			elif (opt_problem._variables[key].type == 'i'):
				raise IOError('NLPQLP cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('NLPQLP cannot handle discrete design variables')
			#end
			i += 1
		#end
		
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
		gg = numpy.zeros([ncon,nproc], numpy.float)
		if ncon > 0:
			i = 0
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					neqc += 1
				#end
				if (self.pll == False):
					gg[i,0] = opt_problem._constraints[key].value
				else:
					for proc in xrange(nproc):
						gg[i,proc] = opt_problem._constraints[key].value
					#end
				#end
				i += 1
			#end
		#end
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		if (len(opt_problem._objectives.keys())>1):
			raise IOError('NLPQLP cannot handle multi-objective problems')
		#end
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff*nproc,numpy.float)
		
		
		# Setup argument list values
		ll = numpy.array([nproc], numpy.int)
		mm = numpy.array([ncon], numpy.int)
		me = numpy.array([neqc], numpy.int)
		mmax = numpy.array([max(1,mm)], numpy.int)
		nn = numpy.array([nvar], numpy.int)
		nmax = numpy.array([max(2,nn+1)], numpy.int)
		mnn2 = numpy.array([mmax+nmax+nmax+2], numpy.int)
		gg = numpy.zeros([mmax,ll], numpy.float)
		df = numpy.zeros([nmax], numpy.float)
		dg = numpy.zeros([mmax,nmax], numpy.float)
		uu = numpy.zeros([mnn2], numpy.float)
		cc = numpy.zeros([nmax,nmax], numpy.float)
		dd = numpy.zeros([nmax], numpy.float)
		acc = numpy.array([self.options['ACC'][1]], numpy.float)
		accqp = numpy.array([self.options['ACCQP'][1]], numpy.float)
		stpmin = numpy.array([self.options['STPMIN'][1]], numpy.float)
		maxfun = numpy.array([self.options['MAXFUN'][1]], numpy.int)
		maxit = numpy.array([self.options['MAXIT'][1]], numpy.int)
		if (nproc != 1):
			maxnm = numpy.array([max(1,min(50,nproc))], numpy.int)
		else:
			maxnm = numpy.array([0], numpy.int)
		#end
		rhob = self.options['RHOB'][1]
		if (self.options['MODE'][1]>=0 and self.options['MODE'][1]<=18):
			mode = self.options['MODE'][1]
		else:
			raise IOError('Incorrect Mode Setting')
		#end
		ifail = numpy.array([0], numpy.int)
		if (myrank == 0):
			if (self.options['IPRINT'][1]>=0 and self.options['IPRINT'][1]<=4):
				iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
			else:
				raise IOError('Incorrect Output Level Setting')
			#end
		else:
			iprint = numpy.array([0], numpy.int)
		#end
		iout = self.options['IOUT'][1]
		ifile = self.options['IFILE'][1]
		if (iprint > 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		lwa = numpy.array([23*nmax+4*mmax+3*mmax+150+3*nmax*nmax/2+10*nmax+nmax+mmax+1], numpy.int)
		wa = numpy.zeros([lwa], numpy.float)
		lkwa = numpy.array([25+nmax], numpy.int)
		kwa = numpy.zeros([lkwa], numpy.intc)
		lactiv = numpy.array([2*mmax+10], numpy.int)
		active = numpy.zeros([lactiv], numpy.bool)
		lql = numpy.array([self.options['LQL'][1]], numpy.bool)
		nfun = numpy.array([0], numpy.int)
		ngrd = numpy.array([0], numpy.int)
		
		# Run NLPQLP
		t0 = time.time()
		xx,ff,gg,uu,nfun,ngrd = nlpqlp.nlpqlp_wrap(ll,mm,me,mmax,nn,nmax,mnn2,xx,ff,gg,df,dg,uu,xl,xu,cc,dd,acc,accqp,stpmin,maxfun,maxit,maxnm,rhob,mode,ifail,iprint,iout,ifile,wa,lwa,kwa,lkwa,active,lactiv,lql,nfun,ngrd,nlfunc,nlgrad)
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
			nlpqlp.closeunit(self.options['IOUT'][1])
		#end
		
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = ifail[0]
		sol_inform['text'] = self.getInform(ifail[0])
		
		if store_sol:
			
			sol_name = 'NLPQLP Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = nfun + ngrd*nvar
			
			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = xx[i,0]
				i += 1
			#end
			
			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = ff[0]
				i += 1
			#end
			
			if ncon > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = -gg[i,0]
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
		
		# 
		if (infocode <= 10):
			return self.informs[infocode]
		else:
			return self.informs[100]
		#end
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  November. 30, 2010 - Ruben E. Perez
		'''
		
		# 
		iprint = self.options['IPRINT'][1]
		if (iprint > 0):
			nlpqlp.pyflush(self.options['IOUT'][1])
		#end
	


#==============================================================================
# NLPQLP Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test NLPQLP
	print 'Testing ...'
	nlpqlp = NLPQLP()
	print nlpqlp
	
