'''
pyMMFD - A Python pyOpt interface to MMFD (part of NASA's ADS).

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.2   $Date: 21/06/2010 21:00$


Tested on:
---------
Win32 with gfortran
Linux with pathf95

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2010)
	v. 1.1	- History support (PJ,RP, 2010)
	v. 1.2  - Gradient Class Support (PJ,RP, 2010)
'''

__version__ = '$Revision: $'

'''
To Do:
	- Include IOUT
	- Implement Informs
	- add unconstrained problems support
'''

# =============================================================================
# MMFD Library
# =============================================================================
try:
	from . import mmfd
except:
	raise ImportError('MMFD shared library failed to import')

import copy
# =============================================================================
# Standard Python modules
# =============================================================================
import os
import sys
import time

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Gradient, Optimizer

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity
# =============================================================================
eps = 1.0	# define a value for machine precision
while ((eps/2.0 + 1.0) > 1.0):
	eps = eps/2.0
eps = 2.0*eps
#eps = math.ldexp(1,-52)


# =============================================================================
# MMFD Optimizer Class
# =============================================================================
class MMFD(Optimizer):

	'''
	MMFD Optimizer Class - Inherited from Optimizer Abstract Class
	'''

	def __init__(self, pll_type=None, *args, **kwargs):

		"""MMFD Optimizer Class Initialization.

		**Keyword arguments:**

		- pll_type -> STR: Parallel Implementation (None, 'POA'-Parallel Objective Analysis), *Default* = None

		Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
		"""

		#
		if (pll_type == None):
			self.poa = False
		elif (pll_type.upper() == 'POA'):
			self.poa = True
		else:
			raise ValueError("pll_type must be either None or 'POA'")

		#
		name = 'MMFD'
		category = 'Local Optimizer'
		def_opts = {
		'IOPT':[int,0],           	# Feasible Directions Approach (0 - MMFD, 1 - MFD)
		'IONED':[int,0],          	# One-Dimensional Search Method (0,1,2,3)
		'CT':[float,-3e-2],       	# Constraint Tolerance
		'CTMIN':[float,4e-3],     	# Active Constraint Tolerance
		'DABOBJ':[float,1e-3],    	# Objective Absolute Tolerance (DABOBJ*abs(f(x)))
		'DELOBJ':[float,1e-3],    	# Objective Relative Tolerance
		'THETAZ':[float,1e-1],    	# Push-Off Factor
		'PMLT':[float,1e1],       	# Penalty multiplier for equality constraints
		'ITMAX':[int,4e2],        	# Maximum Number of Iterations
		'ITRMOP':[int,3],         	# consecutive Iterations Iterations for Convergence
		'IPRINT':[int,2],         	# Print Control (0 - None, 1 - Final, 2 - Iters)
		'IFILE':[str,'MMFD.out'],	# Output File Name
		}
		informs = {
		}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)


	def __solve__(self, opt_problem={}, sens_type='FD', store_sol=True, store_hst=False, hot_start=False, disp_opts=False, sens_mode='', sens_step={}, *args, **kwargs):

		"""Run Optimizer (Optimize Routine)

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
		"""

		#
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pyMMFD - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")

		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print('pyMMFD: Parallel objective Function Analysis requires mpi4py')
			comm = MPI.COMM_WORLD
			nproc = comm.Get_size()
			if (mpi4py.__version__[0] == '0'):
				Bcast = comm.Bcast
			elif (mpi4py.__version__[0] == '1'):
				Bcast = comm.bcast
			self.pll = True
			self.myrank = comm.Get_rank()
		else:
			self.pll = False
			self.myrank = 0

		myrank = self.myrank

		#
		def_fname = self.options['IFILE'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)

		#
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)


		#======================================================================
		# MMFD - Objective/Constraint Values Function
		#======================================================================
		def mmfdfun(nv,nc,x,f,g):

			# Variables Groups Handling
			if opt_problem.use_groups:
				xg = {}
				for group in group_ids.keys():
					if (group_ids[group][1]-group_ids[group][0] == 1):
						xg[group] = x[group_ids[group][0]]
					else:
						xg[group] = x[group_ids[group][0]:group_ids[group][1]]
				xn = xg
			else:
				xn = x

			# Flush Output Files
			self.flushFiles()

			# Evaluate User Function
			fail = 0
			ff = []
			gg = []
			if (myrank == 0):
				if self.hot_start:
					[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
					if hist_end:
						self.hot_start = False
						hos_file.close()
					else:
						[ff,gg,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]

			if self.pll:
				self.hot_start = Bcast(self.hot_start,root=0)
			if self.hot_start and self.pll:
				[ff,gg,fail] = Bcast([ff,gg,fail],root=0)
			elif not self.hot_start:
				[ff,gg,fail] = opt_problem.obj_fun(xn, *args, **kwargs)

			# Store History
			if (myrank == 0):
				if self.sto_hst:
					log_file.write(x,'x')
					log_file.write(ff,'obj')
					log_file.write(gg,'con')
					log_file.write(fail,'fail')

			# Objective Assigment
			if isinstance(ff,complex):
				f = ff.astype(float)
			else:
				f = ff

			# Constraints Assigment
			for i in range(len(opt_problem._constraints.keys())):
				if isinstance(gg[i],complex):
					g[i] = gg[i].astype(float)
				else:
					g[i] = gg[i]

			return f,g


		#======================================================================
		# MMFD - Objective/Constraint Gradients Function
		#======================================================================
		def mmfdgrd(nv,nc,x,f,g,df,dg):

			if self.hot_start:
				dff = []
				dgg = []
				if (myrank == 0):
					[vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
					if hist_end:
						self.hot_start = False
						hos_file.close()
					else:
						dff = vals['grad_obj'][0].reshape((len(opt_problem._objectives.keys()),len(opt_problem._variables.keys())))
						dgg = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))
				if self.pll:
					self.hot_start = Bcast(self.hot_start,root=0)
				if self.hot_start and self.pll:
					[dff,dgg] = Bcast([dff,dgg],root=0)

			if not self.hot_start:

				#
				dff,dgg = gradient.getGrad(x, group_ids, [f], g, *args, **kwargs)


			# Store History
			if self.sto_hst and (myrank == 0):
				log_file.write(dff,'grad_obj')
				log_file.write(dgg,'grad_con')

			# Gradient Assignment
			for i in range(len(opt_problem._variables.keys())):
				df[i] = dff[0,i]
				for j in range(len(opt_problem._constraints.keys())):
					dg[i,j] = dgg[j,i]

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
				raise OSError('MMFD cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise OSError('MMFD cannot handle discrete design variables')
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

		# Constraints Handling
		ncon = len(opt_problem._constraints.keys())
		neqc = 0
		gg = []
		idg = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'i':
					idg.append(1)
				elif opt_problem._constraints[key].type == 'e':
					idg.append(-1)
				gg.append(opt_problem._constraints[key].value)
		else:
			raise OSError('MMFD support for unconstrained problems not implemented yet')
		gg = numpy.array(gg)
		idg = numpy.array(idg, int)

		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		ff = numpy.array(ff)


		# Setup argument list values
		ndv = numpy.array([nvar], int)
		ncn = numpy.array([ncon], int)
		if (self.options['IOPT'][1]>=0 and self.options['IOPT'][1]<=1):
			iopt = numpy.array([self.options['IOPT'][1]], int)
		else:
			raise OSError('Incorrect Feasible Directions Approach')
		if (self.options['IONED'][1]>=0 and self.options['IONED'][1]<=3):
			ioned = numpy.array([self.options['IONED'][1]], int)
		else:
			raise OSError('Incorrect One-Dimensional Search Method')
		if (myrank == 0):
			if (self.options['IPRINT'][1]>=0 and self.options['IPRINT'][1]<=2):
				iprint = numpy.array([self.options['IPRINT'][1]], int)
			else:
				raise OSError('Incorrect Output Level Setting')
		else:
			iprint = numpy.array([0], int)
		#iout = numpy.array([self.options['IOUT'][1]], int)
		ifile = self.options['IFILE'][1]
		if (iprint > 0):
			if os.path.isfile(ifile):
				os.remove(ifile)
		ct = numpy.array([self.options['CT'][1]], float)
		ctmin = numpy.array([self.options['CTMIN'][1]], float)

		finit,ginit = mmfdfun([],[],xx,ff,gg)
		dabobj = numpy.array([self.options['DABOBJ'][1]*finit], float)

		delobj = numpy.array([self.options['DELOBJ'][1]], float)
		thetaz = numpy.array([self.options['THETAZ'][1]], float)
		pmlt = numpy.array([self.options['PMLT'][1]], float)
		itmax = numpy.array([self.options['ITMAX'][1]], int)
		itrmop = numpy.array([self.options['ITRMOP'][1]], int)
		nrwk0 = 500
		nrwk1 = 10*(2*nvar+ncon)
		nrwk2 = (ncon+2*nvar+3)
		nrwk3 = (ncon+2*nvar)*((ncon+2*nvar)/2+1)
		nrwkS = nrwk0 + nrwk1 + nrwk2 + nrwk3
		nrwk = numpy.array([nrwkS], int)
		wk = numpy.zeros([nrwk], float)
		nriwk = numpy.array([nrwkS], int)
		iwk = numpy.zeros([nriwk], int)

		nfun = numpy.array([0], int)
		ngrd = numpy.array([0], int)


		# Run MMFD
		t0 = time.time()
		mmfd.mmfd(iopt,ioned,iprint,ndv,ncn,xx,xl,xu,ff,gg,idg,
			wk,nrwk,iwk,nriwk,ifile,ct,ctmin,dabobj,delobj,thetaz,
			pmlt,itmax,itrmop,nfun,ngrd,mmfdfun,mmfdgrd)
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

		if (iprint > 0):
		#	mmfd.closeunit(self.options['IOUT'][1])
			mmfd.closeunit(6)


		# Store Results
		sol_inform = {}
		sol_inform['value'] = []
		sol_inform['text'] = {}

		if store_sol:

			sol_name = 'MMFD Solution to ' + opt_problem.name

			sol_options = copy.copy(self.options)
			if 'defaults' in sol_options:
				del sol_options['defaults']

			sol_evals = nfun[0] + ngrd[0]*nvar

			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = xx[i]
				i += 1

			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = ff[i]
				i += 1

			if ncon > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = gg[i]
					i += 1
			else:
				sol_cons = {}

			sol_lambda = {}


			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time,
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options,
				display_opts=disp_opts, Lambda=sol_lambda, Sensitivities=sens_type,
				myrank=myrank, arguments=args, **kwargs)


		return ff, xx, sol_inform



	def _on_setOption(self, name, value):

		"""Set Optimizer Option Value (Optimizer Specific Routine)

		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		"""

		pass


	def _on_getOption(self, name):

		"""Get Optimizer Option Value (Optimizer Specific Routine)

		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		"""

		pass


	def _on_getInform(self, infocode):

		"""Get Optimizer Result Information (Optimizer Specific Routine)

		Keyword arguments:
		-----------------
		id -> STRING: Option Name

		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		"""

		pass


	def _on_flushFiles(self):

		"""Flush the Output Files (Optimizer Specific Routine)

		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		"""

		#
		iPrint = self.options['IPRINT'][1]
		if (iPrint > 0):
			#mmfd.pyflush(self.options['IOUT'][1])
			mmfd.pyflush(6)
