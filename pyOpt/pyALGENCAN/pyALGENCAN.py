#!/usr/bin/env python
'''
pyALGENCAN - A Python pyOpt interface to ALGENCAN. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.0   $Date: 10/06/2014 21:00$


Tested on:
---------
Win32 with g77
Linux with g77

Developers:
-----------
- Dr. Ruben E. Perez (RP)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2011)
	V. 1.1  - Updated to Work with ALGENCAN 2.4 (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- add evals return on source
	- add unconstrained problems support
'''

# =============================================================================
# ALGENCAN Library
# =============================================================================
try:
	import algencan
except ImportError:
	raise ImportError('ALGENCAN shared library failed to import')
#end

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import copy, time
import pdb

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
# ALGENCAN Optimizer Class
# =============================================================================
class ALGENCAN(Optimizer):
	
	'''
	ALGENCAN Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		ALGENCAN Optimizer Class Initialization
		
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
		name = 'ALGENCAN'
		category = 'Local Optimizer'
		def_opts = {
		# ALGENCAN Options
		'epsfeas':[float,1.0e-8],		# Feasibility Convergence Accurancy
		'epsopt':[float,1.0e-8],		# Optimality Convergence Accurancy
		'efacc':[float,1.0e-4],			# Feasibility Level for Newton-KKT Acceleration
		'eoacc':[float,1.0e-4],			# Optimality Level for Newton-KKT Acceleration
		'checkder':[bool,False],		# Check Derivatives Flag
		'iprint':[int,10],				# Print Flag (0 - None, )
		'ifile':[str,'ALGENCAN.out'],	# Output File Name
		'ncomp':[int,6],				# Print Precision
		}
		informs = {
		0 : "Solution was found.",
		1 : "Stationary or infeasible point was found.",
		2 : "penalty parameter is too large infeasibile or badly scaled problem",
		3 : "Maximum of iterations reached.",
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
			raise NotImplementedError("pyALGENCAN - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyALGENCAN: Parallel objective Function Analysis requires mpi4py'
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
		
		
		# 
		def evalf(n,x,f,flag):
			flag = -1
			return f,flag
		
		def evalg(n,x,g,flag):
			flag = -1
			return g,flag
		
		def evalh(n,x,hlin,hcol,hval,hnnz,flag):
			flag = -1
			return hlin,hcol,hval,hnnz,flag
		
		def evalc(n,x,ind,c,flag):
			flag = -1
			return c,flag
		
		def evaljac(n,x,ind,jcvar,jcval,jcnnz,flag):
			flag = -1
			return jcvar,jcval,jcnnz,flag
		
		def evalgjacp(n,x,g,m,p,q,work,gotj,flag):
			flag = -1
			return g,p,q,gotj,flag
		
		def evalhc(n,x,ind,hclin,hccol,hcval,hcnnz,flag):
			flag = -1
			return hclin,hccol,hcval,hcnnz,flag
		
		def evalhl(n,x,m,lmbda,scalef,scalec,hllin,hlcol,hlval,hlnnz,flag):
			flag = -1
			return hllin,hlcol,hlval,hlnnz,flag
		
		def evalhlp(n,x,m,lmbda,sf,sc,p,hp,goth,flag):
			flag = -1
			return hp,goth,flag
		
		
		#======================================================================
		# ALGENCAN - Objective/Constraint Values Function
		#======================================================================
		def evalfc(n,x,f,m,g,flag):
			
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
			flag = 0
			if (myrank == 0):
				if self.h_start:
					[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						[ff,gg,flag] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
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
			
			# Constraints Assigment
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(gg[i],complex):
					g[i] = gg[i].astype(float)
				else:
					g[i] = gg[i]
				#end
			#end
			
			return f,g,fail
		
		
		#======================================================================
		# ALGENCAN - Objective/Constraint Gradients Function
		#======================================================================
		def evalgjac(n,x,jfval,m,jcfun,jcvar,jcval,jcnnz,flag):
			
			if self.h_start:
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
				
				[ff,gg,fail] = opt_problem.obj_fun(x, *args, **kwargs)
				dff,dgg = gradient.getGrad(x, group_ids, [ff], gg, *args, **kwargs)
				
			#end
			
			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')
			#end
			
			# Objective Gradient Assigment
			for i in xrange(len(opt_problem._variables.keys())):
				jfval[i] = dff[0,i]
			#end
			
			# Constraint Gradient Assigment
			jcnnz = 0
			for jj in xrange(len(opt_problem._constraints.keys())):
				for ii in xrange(len(opt_problem._variables.keys())):
					jcfun[jcnnz] = jj + 1
					jcvar[jcnnz] = ii + 1
					jcval[jcnnz] = dgg[jj,ii]
					jcnnz += 1
				#end
			#end
			
			return jfval,jcfun,jcvar,jcval,jcnnz,fail
		
		
		
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
		m = len(opt_problem._constraints.keys())
		equatn = []
		linear = []
		if m > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					equatn.append(True)
				elif opt_problem._constraints[key].type == 'i':
					equatn.append(False)
				#end
				linear.append(False)
			#end
		else:
			raise IOError('ALGENCAN support for unconstrained problems not implemented yet')
		#end
		equatn = numpy.array(equatn)
		linear = numpy.array(linear)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff)
		
		
		# Setup argument list values
		nn = numpy.array([n], numpy.int)
		mm = numpy.array([m], numpy.int)
		lm = numpy.zeros([m], numpy.float)
		coded = numpy.array([False,False,False,False,False,False,True,True,False,False], numpy.bool)
		epsfeas = numpy.array([self.options['epsfeas'][1]], numpy.float)
		epsopt = numpy.array([self.options['epsopt'][1]], numpy.float)
		efacc = numpy.array([self.options['efacc'][1]], numpy.float)
		eoacc = numpy.array([self.options['eoacc'][1]], numpy.float)	
		checkder = numpy.array([self.options['checkder'][1]], numpy.bool)
		iprint = numpy.array([self.options['iprint'][1]], numpy.int)
		if (myrank != 0):
			iprint = 0
		else:
			iprint = self.options['iprint'][1]
		#end
		ncomp = numpy.array([self.options['ncomp'][1]], numpy.int)

		ifile = self.options['ifile'][1]
		if (iprint >= 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		cnormu = numpy.array([0], numpy.float)
		snorm = numpy.array([0], numpy.float)
		nlpsupn = numpy.array([0], numpy.float)
		inform = numpy.array([0], numpy.int)
		
		# Run ALGENCAN
		t0 = time.time()
		algencan.algencan(epsfeas,epsopt,efacc,eoacc,iprint,ncomp,nn,xx,xl,xu,mm,lm,equatn,linear,coded,checkder,ff,cnormu,snorm,nlpsupn,inform,ifile,evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalgjacp,evalhl,evalhlp)
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
			algencan.closeunit(10)
		#end
		
		
		# 
		[fs,gg,fail] = opt_problem.obj_fun(xx, *args, **kwargs)
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = inform[0]
		sol_inform['text'] = self.getInform(inform[0])
		
		if store_sol:
			
			sol_name = 'ALGENCAN Solution to ' + opt_problem.name
			
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
					sol_cons[key].value = gg[i]
					i += 1
				#end
			else:
				sol_cons = {}
			#end
			
			sol_lambda = lm
			
			
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
		iPrint = self.options['iprint'][1]
		if (iPrint >= 0):
			algencan.pyflush(10)	
		#end
	


#==============================================================================
# ALGENCAN Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test ALGENCAN
	print 'Testing ...'
	algencan = ALGENCAN()
	print algencan
	
