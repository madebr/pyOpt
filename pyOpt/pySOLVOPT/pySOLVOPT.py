#!/usr/bin/env python
'''
pySOLVOPT - A Python pyOpt interface to SOLVOPT. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.5   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with g77
Linux with pathf95
Win32 with g77

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Andrew Lambe (AL)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0  - Initial Class Creation (RP, 2009)
	v. 1.1	- Integrate to pyOpt Framework (RP, 2009)
	v. 1.2  - Wrapper Callback Storage Support (AL,RP, 2009)
	v. 1.3  - History support (PJ,RP, 2010)
	v. 1.4  - Gradient Class Support (PJ,RP, 2010)
	v. 1.5  - Unconstrained Problems Support (RP, 2014) 
'''

__version__ = '$Revision: $'

'''
To Do:
	- check possible ff,gg issue in eval
'''

# =============================================================================
# SOLVOPT Library
# =============================================================================
try:
	import solvopt
except:
	raise ImportError('SOLVOPT shared library failed to import')
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
# SOLVOPT Optimizer Class
# =============================================================================
class SOLVOPT(Optimizer):
	
	'''
	SOLVOPT Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		SOLVOPT Optimizer Class Initialization
		
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
		name = 'SOLVOPT'
		category = 'Local Optimizer'
		def_opts = {
		'xtol':[float,1e-4],			# Variables Tolerance
		'ftol':[float,1e-6],			# Objective Tolerance
		'maxit':[int,15000],			# Maximum Number of Iterations
		'iprint':[int,1],     			# Output Level (-1 -> None, 0 -> Final, N - each Nth iter)
		'gtol':[float,1e-8],  			# Constraints Tolerance
		'spcdil':[float,2.5], 			# Space Dilation 
		'iout':[int,6],     			# Output Unit Number
		'ifile':[str,'SOLVOPT.out'],	# Output File Name
		}
		informs = {
		1 : 'Normal termination.', 
		-2 : 'Improper space dimension.',
		-3 : 'Objective equals infinity.',
		-4 : 'Gradient equals zero or infinity.',
		-5 : 'Objective equals infinity.',
		-6 : 'Gradient equals zero or infinity.',
		-7 : 'Objective function is unbounded.',
		-8 : 'Gradient zero at the point, but stopping criteria are not fulfilled.',
		-9 : 'Iterations limit exceeded.',
		-11 : 'Premature stop is possible. Try to re-run the routine from the obtained point.',
		-12 : 'Result may not provide the optimum. The function apparently has many extremum points.',
		-13 : 'Result may be inaccurate in the coordinates. The function is flat at the optimum.',
		-14 : 'Result may be inaccurate in a function value. The function is extremely steep at the optimum.',
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
		
		Documentation last updated:  February. 2, 2011 - Peter W. Jansen
		'''
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pySOLVOPT - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pySOLVOPT: Parallel objective Function Analysis requires mpi4py'
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
		# SOLVOPT - Objective/Constraint Values Storage 
		#======================================================================
		def soeval(x):
			
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
				df,dg = gradient.getGrad(x, group_ids, [f], g, *args, **kwargs)
				
			#end
			
			# Store History
			if (myrank == 0):
				if self.sto_hst:
					log_file.write(df,'grad_obj')
					log_file.write(dg,'grad_con')
				#end
			#end
			
			# Objective Assigment
			if isinstance(f,complex):
				ff = f.astype(float)
			else:
				ff = f
			#end
			
			# Constraints Assigment
			i = 0
			j = 0
			for j in xrange(len(opt_problem._constraints.keys())):
				if isinstance(g[j],complex):
					gg[i] = g[j].astype(float)
				else:
					gg[i] = g[j]
				#end
				i += 1
			#end
			for key in opt_problem._variables.keys():
				if (opt_problem._variables[key].lower != -inf):
					gg[i] = xl[j] - x[j]
					i += 1
				#end
				if (opt_problem._variables[key].upper != inf):
					gg[i] = x[j] - xu[j]
					i += 1
				#end
			#end
			
			# Gradient Assignment
			df = df[0]
			dgg = numpy.zeros([len(opt_problem._variables.keys())*2,len(opt_problem._variables.keys())],'d')
			i = len(opt_problem._constraints.keys())
			j = 0
			for key in opt_problem._variables.keys():
				if (opt_problem._variables[key].lower != -inf): 
					if (gg[i] > 0):
						dgg[j,j] = -1
					#end	
					i += 1
				#end
				if (opt_problem._variables[key].upper != inf): 
					if (gg[i] > 0):
						dgg[j,j] = 1
					#end
					i += 1
				#end
				j += 1
			#end
			dg = numpy.concatenate((dg,dgg),axis=0)
			
			
			# Store
			self.stored_data['x'] = copy.copy(x)
			self.stored_data['f'] = copy.copy(ff)
			self.stored_data['g'] = copy.copy(gg)
			self.stored_data['df'] = copy.copy(df)
			self.stored_data['dg'] = copy.copy(dg)
			
			return
		
		
		#======================================================================
		# SOLVOPT - Objective Value Function
		#======================================================================
		def soobjf(n,x,f):
			
			if ((self.stored_data['x'] != x).any()):
				soeval(x)
			#end
			
			f = self.stored_data['f']
			
			return f
		
		
		#======================================================================
		# SOLVOPT - Constraint Values Function
		#======================================================================
		def soobjg(n,x,g):
			
			if ((self.stored_data['x'] != x).any()):
				soeval(x)
			#end
			
			# Constraints Maximal Residual
			maxg = numpy.zeros([len(self.stored_data['g'])],float)
			i = 0
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					maxg[i] = abs(self.stored_data['g'][i])
					i += 1
				elif opt_problem._constraints[key].type == 'i':
					maxg[i] = max(0,self.stored_data['g'][i])
					i += 1
				#end
			#end
			for j in xrange(len(opt_problem._constraints),len(self.stored_data['g'])):
				maxg[i] = max(0,self.stored_data['g'][i])
				i += 1
			#end
			
			g = max(maxg)
			
			return g
		
		
		#======================================================================
		# SOLVOPT - Objective Gradients Function
		#======================================================================
		def sogrdf(n,x,df):
			
			if ((self.stored_data['x'] != x).any()):
				soeval(x)
			#end
			
			df = self.stored_data['df']
			
			return df
		
		
		#======================================================================
		# SOLVOPT - Constraint Gradients Function
		#======================================================================
		def sogrdg(n,x,dg):
			
			if ((self.stored_data['x'] != x).any()):
				soeval(x)
			#end
			
			# Constraints Maximal Residual
			maxg = numpy.zeros([len(self.stored_data['g'])],float)
			i = 0
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					maxg[i] = abs(self.stored_data['g'][i])
					i += 1
				elif opt_problem._constraints[key].type == 'i':
					maxg[i] = max(0,self.stored_data['g'][i])
					i += 1
				#end
			#end
			for j in xrange(len(opt_problem._constraints),len(self.stored_data['g'])):
				maxg[i] = max(0,self.stored_data['g'][i])
				i += 1
			#end
			id = min(numpy.nonzero(maxg==max(maxg))[0])
			
			dg = self.stored_data['dg'][id]
			
			return dg
		
		
		
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
				raise IOError('SOLVOPT cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('SOLVOPT cannot handle discrete design variables')
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
		gg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					neqc += 1
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
		flg = numpy.array([True], numpy.bool)
		iprint = self.options['iprint'][1]
		if (myrank != 0):
			iprint = -1
		else:
			iprint = self.options['iprint'][1]
		#end
		options = numpy.zeros([13], numpy.float)
		options[0] = -1  						# Minimize
		options[1] = self.options['xtol'][1]	# Variables Tolerance
		options[2] = self.options['ftol'][1]	# Objective Tolerance
		options[3] = self.options['maxit'][1]	# Maximum Number of Iterations
		options[4] = iprint						# Output Level
		options[5] = self.options['gtol'][1]	# Constraints Tolerance
		options[6] = self.options['spcdil'][1]	# Space Dilation
		options[7] = 1e-11						# LB FD Stepsize (NA as we provide our own sensitivities)
		flfc = numpy.array([True], numpy.bool)
		flgc = numpy.array([True], numpy.bool)
		iout = numpy.array([self.options['iout'][1]], numpy.int)
		ifile = self.options['ifile'][1]
		
		if (iprint >= 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
			#end
		#end
		wb = numpy.zeros([nvar,nvar], numpy.float)
		wg = numpy.zeros([nvar], numpy.float)
		wg0 = numpy.zeros([nvar], numpy.float)
		wg1 = numpy.zeros([nvar], numpy.float)
		wgt = numpy.zeros([nvar], numpy.float)
		wgc = numpy.zeros([nvar], numpy.float)
		wz = numpy.zeros([nvar], numpy.float)
		wx1 = numpy.zeros([nvar], numpy.float)
		wxopt = numpy.zeros([nvar], numpy.float)
		wxrec = numpy.zeros([nvar], numpy.float)
		wgrec = numpy.zeros([nvar], numpy.float)
		wxx = numpy.zeros([nvar], numpy.float)
		wdeltax = numpy.zeros([nvar], numpy.float)
		widx = numpy.zeros([nvar], numpy.int)
		
		
		# Storage Arrays
		self.stored_data = {}
		self.stored_data['x'] = {}  #numpy.zeros([nvar],float)
		self.stored_data['f'] = {}  #numpy.zeros([nobj],float)
		self.stored_data['g'] = {}  #numpy.zeros([ncon+nadd],float)
		self.stored_data['df'] = {} #numpy.zeros([nvar],float)
		self.stored_data['dg'] = {} #numpy.zeros([ncon+nadd,nvar],float) 
		
		
		# Run SOLVOPT
		t0 = time.time()
		solvopt.solvopt(n,xx,ff,soobjf,flg,sogrdf,options,flfc,soobjg,flgc,sogrdg,wb,wg,wg0,wg1,wgt,wgc,wz,wx1,wxopt,wxrec,wgrec,wxx,wdeltax,widx,iout,ifile)
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
			solvopt.closeunit(self.options['iout'][1])
		#end
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = options[8]
		sol_inform['text'] = self.getInform(options[8])
		
		if store_sol:
			
			sol_name = 'SOLVOPT Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = options[9]+options[10] # assumes cnst fevals and gevals are included in feval & geval
			
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
		
		if (infocode > 0):
			return self.informs[1]
		else:
			return self.informs[infocode]
		#end
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		'''
		
		# 
		iprint = self.options['iprint'][1]
		if (iprint >= 0):
			solvopt.pyflush(self.options['iout'][1])	
		#end
	


#==============================================================================
# SOLVOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test SOLVOPT
	print 'Testing ...'
	solvopt = SOLVOPT()
	print solvopt
	
