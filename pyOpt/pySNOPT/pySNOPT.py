#! /usr/bin/env python
'''
pySNOPT - A Python pyOpt interface to SNOPT.

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.4   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with g77
Win32 with g77
Linux with intel9
Linux with intel11
Linux with pathf95
Linux with gfortran
Mac with g95
scinet gpc intel 11 (--fcompiler=intelem)

Developers:
-----------
- Dr. Joaquim Martins (JM)
- Dr. Patrick LeGresley (PL)
- Mr. Nathan Tedford (NT)
- Mr. C.A.(Sandy) Mader (SM)
- Mr. Chris Marriage (CM)
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
	v. 0.1	- Initial Wrapper Creation (JM, 2000)
	v. 0.2	- Wrapper Bug Fixes and Updates (PL, 2003-2005)
	v. 0.3	- Wrapper Additional Functionality Updates (NT, 2006)
	v. 0.4	- Wrapper Updates to Work with NumPy (SM, 2007)
	v. 0.5	- Wrapper Updates to Work with SNOPT 7.0 (CM, 2007)
	v. 1.0 	- Initial pyOpt Class Creation (RP, 2008)
			- Migrate Wrapper to pyOpt Framework (RP, 2008)
	v. 1.1  - User provided sensitivities support (PJ,RP, 2008)
	v. 1.2	- History support (PJ,RP, 2010)
	v. 1.3	- Gradient Class Support (PJ,RP, 2010)
	v. 1.4  - Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- check how to properly do internal snopt sensitivities!
	- snopt mode flag handling!
	- number of evals from where?
	- read/write spec to/from options!
	- check x from history for deviations in hot start
'''

# =============================================================================
# SNOPT Library
# =============================================================================
try:
	import snopt
except:
	raise ImportError('SNOPT shared library failed to import')
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
# SNOPT Optimizer Class
# =============================================================================
class SNOPT(Optimizer):
	
	'''
	SNOPT Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		SNOPT Optimizer Class Initialization
		
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
		name = 'SNOPT'
		category = 'Local Optimizer'
		def_opts = {
		# SNOPT Printing Options
		'Major print level':[int,1], 					# Majors Print (1 - line major iteration log)
		'Minor print level':[int,1], 					# Minors Print (1 - line minor iteration log)
		'Print file':[str,'SNOPT_print.out'],			# Print File Name (specified by subroutine snInit)
		'iPrint':[int,18],								# Print File Output Unit (override internally in snopt?)
		'Summary file':[str,'SNOPT_summary.out'],		# Summary File Name (specified by subroutine snInit)
		'iSumm':[int,19],								# Summary File Output Unit (override internally in snopt?)
		'Print frequency':[int,100],					# Minors Log Frequency on Print File
		'Summary frequency':[int,100],					# Minors Log Frequency on Summary File
		'Solution':[str,'Yes'],							# Print Solution on the Print File
		'Suppress options listing':[type(None),None],	# (options are normally listed)
		'System information':[str,'No'],				# Print System Information on the Print File
		# SNOPT Problem Specification Options
		'Problem Type':[str,'Minimize'],				# ('Maximize': alternative over Minimize, 'Feasible point': alternative over Minimize or Maximize)
		'Objective row':[int,1], 						# (has precedence over ObjRow (snOptA))
		'Infinite bound':[float,1.0e+20],				# Infinite Bound Value
		# SNOPT Convergence Tolerances Options
		'Major feasibility tolerance':[float,1.0e-6],	# Target Nonlinear Constraint Violation
		'Major optimality tolerance':[float,1.0e-6], 	# Target Complementarity Gap
		'Minor feasibility tolerance':[float,1.0e-6], 	# For Satisfying the QP Bounds
		# SNOPT Derivative Checking Options
		'Verify level':[int,0],							# Gradients Check Flag
		# SNOPT Scaling Options
		'Scale option':[int,1],							# Scaling (1 - linear constraints and variables)
		'Scale tolerance':[float,0.9],					# Scaling Tolerance
		'Scale Print':[type(None),None],				# Default: scales are not printed
		# SNOPT Other Tolerances Options
		'Crash tolerance':[float,0.1],					# 
		'Linesearch tolerance':[float,0.9], 			# smaller for more accurate search
		'Pivot tolerance':[float,3.7e-11], 				# epsilon^(2/3)
		# SNOPT QP subproblems Options
		'QPSolver':[str,'Cholesky'], 					# Default: Cholesky
		'Crash option':[int,3], 						# (3 - first basis is essentially triangular)
		'Elastic mode':[str,'No'],						# (start with elastic mode until necessary)
		'Elastic weight':[float,1.0e+4], 				# (used only during elastic mode)
		'Iterations limit':[int,10000], 				# (or 20*ncons if that is more)
		'Partial price':[int,1], 						# (10 for large LPs)
		# SNOPT SQP method Options
		'Start':[str,'Cold'], 							# has precedence over argument start, ('Warm': alternative to a cold start)
		'Major iterations limit':[int,1000], 			# or ncons if that is more
		'Minor iterations limit':[int,500], 			# or 3*ncons if that is more
		'Major step limit':[float,2.0],					# 
		'Superbasics limit':[int,None],     			# (n1 + 1, n1 = number of nonlinear variables)
		'Derivative level':[int,3],						# (NOT ALLOWED IN snOptA)
		'Derivative option':[int,1],					# (ONLY FOR snOptA)
		'Derivative linesearch':[type(None),None],		#
		'Nonderivative linesearch':[type(None),None],	#
		'Function precision':[float,3.0e-13], 			# epsilon^0.8 (almost full accuracy)
		'Difference interval':[float,5.5e-7], 			# Function precision^(1/2)
		'Central difference interval':[float,6.7e-5],	# Function precision^(1/3)
		'New superbasics limit':[int,99],				# controls early termination of QPs
		'Objective row':[int,1],						# row number of objective in F(x)
		'Penalty parameter':[float,0.0], 				# initial penalty parameter
		'Proximal point method':[int,1],				# (1 - satisfies linear constraints near x0)
		'Reduced Hessian dimension':[int,2000],			# (or Superbasics limit if that is less)
		'Violation limit':[int,10.0], 					# (unscaled constraint violation limit)
		'Unbounded step size':[float,1.0e+18],			# 
		'Unbounded objective':[float,1.0e+15],			# 
		# SNOPT Hessian approximation Options
		'Hessian full memory':[type(None),None], 		# default if n1 <= 75
		'Hessian limited memory':[type(None),None], 	# default if n1 > 75
		'Hessian frequency':[int,999999], 				# for full Hessian (never reset)
		'Hessian updates':[int,10],	 					# for limited memory Hessian
		'Hessian flush':[int,999999], 					# no flushing
		# SNOPT Frequencies Options
		'Check frequency':[int,60],	 					# test row residuals ||Ax - sk||
		'Expand frequency':[int,10000],					# for anti-cycling procedure
		'Factorization frequency':[int,50],				# 100 for LPs
		'Save frequency':[int,100],	 					# save basis map
		# SNOPT LUSOL Options
		'LU factor tolerance':[float,3.99], 			# for NP (100.0 for LP)
		'LU update tolerance':[float,3.99], 			# for NP ( 10.0 for LP)
		'LU singularity tolerance':[float,3.2e-11],		# 
		'LU partial pivoting':[type(None),None], 		# default threshold pivoting strategy
		'LU rook pivoting':[type(None),None], 			# threshold rook pivoting
		'LU complete pivoting':[type(None),None], 		# threshold complete pivoting
		# SNOPT Basis files Options
		'Old basis file':[int,0], 						# input basis map
		'New basis file':[int,0], 						# output basis map
		'Backup basis file':[int,0], 					# output extra basis map
		'Insert file':[int,0], 							# input in industry format
		'Punch file':[int,0], 							# output Insert data
		'Load file':[int,0], 							# input names and values
		'Dump file':[int,0], 							# output Load data
		'Solution file':[int,0], 						# different from printed solution
		# SNOPT Partitions of cw, iw, rw Options
		'Total character workspace':[int,500],  		# lencw: 500
		'Total integer workspace':[int,None],			# leniw: 500 + 100 * (m+n) 
		'Total real workspace':[int,None],				# lenrw: 500 + 200 * (m+n)
		'User character workspace':[int,500],			# 
		'User integer workspace':[int,500],				# 
		'User real workspace':[int,500],				# 
		#SNOPT Miscellaneous Options
		'Debug level':[int,0], 							# (0 - Normal, 1 - for developers)
		'Timing level':[int,3],							# (3 - print cpu times)
		}
		informs = {
		0 : 'finished successfully',
		1 : 'optimality conditions satisfied',
		2 : 'feasible point found',
		3 : 'requested accuracy could not be achieved',
		4 : 'weak QP minimizer',
		10 : 'the problem appears to be infeasible',
		11 : 'infeasible linear constraints',
		12 : 'infeasible linear equalities',
		13 : 'nonlinear infeasibilities minimized',
		14 : 'infeasibilities minimized',
		15 : 'infeasible linear constraints in QP subproblem',
		20 : 'the problem appears to be unbounded',
		21 : 'unbounded objective',
		22 : 'constraint violation limit reached',
		30 : 'resource limit error',
		31 : 'iteration limit reached',
		32 : 'major iteration limit reached',
		33 : 'the superbasics limit is too small',
		40 : 'terminated after numerical difficulties',
		41 : 'current point cannot be improved',
		42 : 'singular basis',
		43 : 'cannot satisfy the general constraints',
		44 : 'ill-conditioned null-space basis',
		50 : 'error in the user-supplied functions',
		51 : 'incorrect objective  derivatives',
		52 : 'incorrect constraint derivatives',
		53 : 'the QP Hessian is indefinite',
		54 : 'incorrect second derivatives',
		55 : 'incorrect derivatives',
		60 : 'undefined user-supplied functions',
		61 : 'undefined function at the first feasible point',
		62 : 'undefined function at the initial point',
		63 : 'unable to proceed into undefined region',
		70 : 'user requested termination',
		71 : 'terminated during function evaluation',
		72 : 'terminated during constraint evaluation',
		73 : 'terminated during objective evaluation',
		74 : 'terminated from monitor routine',
		80 : 'insufficient storage allocated',
		81 : 'work arrays must have at least 500 elements',
		82 : 'not enough character storage',
		83 : 'not enough integer storage',
		84 : 'not enough real storage',
		90 : 'input arguments out of range',
		91 : 'invalid input argument',
		92 : 'basis file dimensions do not match this problem',
		93 : 'the QP Hessian is indefinite',
		100 : 'finished successfully',
		101 : 'SPECS file read',
		102 : 'Jacobian structure estimated',
		103 : 'MPS file read',
		104 : 'memory requirements estimated',
		105 : 'user-supplied derivatives appear to be correct',
		106 : 'no derivatives were checked',
		107 : 'some SPECS keywords were not recognized',
		110 : 'errors while processing MPS data',
		111 : 'no MPS file specified',
		112 : 'problem-size estimates too small',
		113 : 'fatal error in the MPS file',
		120 : 'errors while estimating Jacobian structure',
		121 : 'cannot find Jacobian structure at given point',
		130 : 'fatal errors while reading the SP',
		131 : 'no SPECS file (iSpecs le 0 or iSpecs gt 99)',
		132 : 'End-of-file while looking for a BEGIN',
		133 : 'End-of-file while reading SPECS file',
		134 : 'ENDRUN found before any valid SPECS',
		140 : 'system error',
		141 : 'wrong no of basic variables',
		142 : 'error in basis package',
		142 : 'Problem dimensions are too large'
		}
		self.set_options = []
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
		
		Documentation last updated:  Feb. 2, 2011 - Peter W. Jansen
		'''
		
		# 
		if ((self.poa) and (sens_mode.lower() == 'pgc')):
			raise NotImplementedError("pySNOPT - Current implementation only allows single level parallelization, either 'POA' or 'pgc'")
		#end
		
		if self.poa or (sens_mode.lower() == 'pgc'):
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pySNOPT: Parallel objective Function Analysis requires mpi4py'
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
		def_fname = self.options['Print file'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		#
		gradient = Gradient(opt_problem, sens_type, sens_mode, sens_step, *args, **kwargs)
		
		
		#======================================================================
		# SNOPT-C - Objective/Constraint Values Function
		#======================================================================
		def snfuncgrag(mode, njac, x, f_obj, g_obj, f_con, g_con):
			
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
			if (myrank == 0):
				if self.h_start:
					[vals,hist_end] = hos_file.read(ident=['obj', 'con', 'fail'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						[f_obj,f_con,fail] = [vals['obj'][0][0],vals['con'][0],int(vals['fail'][0][0])]
					#end
				#end
			#end
			
			if self.pll:
				self.h_start = Bcast(self.h_start,root=0)
			#end
			if self.h_start and self.pll:
				[f_obj,f_con,fail] = Bcast([f_obj,f_con,fail],root=0)
			elif not self.h_start:	
				[f_obj,f_con,fail] = opt_problem.obj_fun(xn, *args, **kwargs)
			#end
			
			# Store History
			if (myrank == 0):
				if self.sto_hst:
					log_file.write(x,'x')
					log_file.write(f_obj,'obj')
					log_file.write(f_con,'con')
					log_file.write(fail,'fail')
				#end
			#end
			
			# 
			if (fail == 1):
				mode = -1
				return mode
			#end
			
			# Gradients
			if mode <> 0 and self.h_start:
				if (myrank == 0):
					[vals,hist_end] = hos_file.read(ident=['grad_obj','grad_con'])
					if hist_end:
						self.h_start = False
						hos_file.close()
					else:
						g_obj = vals['grad_obj'][0]
						g_con = vals['grad_con'][0].reshape((len(opt_problem._constraints.keys()),len(opt_problem._variables.keys())))	
					#end
				#end
				if self.pll:
					self.h_start = Bcast(self.h_start,root=0)
				#end
				if self.h_start and self.pll:
					[g_obj,g_con] = Bcast([g_obj,g_con],root=0)
				#end
			#end
			
			if mode <> 0 and not self.h_start:
				
				# 
				dff,dgg = gradient.getGrad(x, group_ids, [f_obj], f_con, *args, **kwargs)
				
				# 
				for i in xrange(len(opt_problem._variables.keys())):
					g_obj[i] = dff[0,i]
					for j in xrange(len(opt_problem._constraints.keys())):
						g_con[j,i] = dgg[j,i]
					#end
				#end
			#end
			
			if (myrank == 0):
				if mode <> 0 and self.sto_hst:
					log_file.write(g_obj,'grad_obj')
					log_file.write(g_con,'grad_con')
				#end
			#end
			
			# Objective Assigment
			if isinstance(f_obj,complex):
				f_obj = f_obj.astype(float)
			#end
			
			# Constraints Assigment
			for i in xrange(len(opt_problem._constraints.keys())):
				if isinstance(f_con[i],complex):
					f_con[i] = f_con[i].astype(float)
				#end
			#end
			if not f_con:
				f_con = [0]
			#end
			
			return mode,f_obj,g_obj,f_con,g_con
		
		
		
		# Variables Handling
		nvar = len(opt_problem._variables.keys())
		blx = []
		bux = []
		xs = []
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'c'):
				blx.append(opt_problem._variables[key].lower)
				bux.append(opt_problem._variables[key].upper)
				xs.append(opt_problem._variables[key].value)
			elif (opt_problem._variables[key].type == 'i'):
				raise IOError('SNOPT cannot handle integer design variables')
			elif (opt_problem._variables[key].type == 'd'):
				raise IOError('SNOPT cannot handle discrete design variables')
			#end
		#end
		blx = numpy.array(blx)
		bux = numpy.array(bux)
		xs = numpy.array(xs)
		
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
		blc = []
		buc = []
		if ncon > 0:
			for key in opt_problem._constraints.keys():
				if (opt_problem._constraints[key].type == 'e'):
					blc.append(opt_problem._constraints[key].equal)
					buc.append(opt_problem._constraints[key].equal)
				elif (opt_problem._constraints[key].type == 'i'):
					blc.append(opt_problem._constraints[key].lower)
					buc.append(opt_problem._constraints[key].upper)
				#end
			#end
		else:
			#if ((store_sol) and (myrank == 0)):
			#	print "Optimization Problem Does Not Have Constraints\n"
			#	print "Unconstrained Optimization Initiated\n"
			##end
			ncon = 1
			blc.append(-inf)
			buc.append( inf)
		#end
		blc = numpy.array(blc)
		buc = numpy.array(buc)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff, numpy.float)
		
		
		# Initialize SNOPT
		iPrint = self.options['iPrint'][1]
		if (myrank != 0):
			iPrint = 0
		#end
		PrintFile = self.options['Print file'][1]
		if (iPrint != 0):
			if os.path.isfile(PrintFile):
				os.remove(PrintFile)
			#end
			ierror = snopt.openunit(iPrint, numpy.array(PrintFile), numpy.array('new'), numpy.array('sequential'))
			if (ierror != 0):
				raise IOError('Failed to properly open %s, ierror = %3d' %(PrintFile,ierror))
			#end
		#end
		
		iSumm = self.options['iSumm'][1]
		if (myrank != 0):
			iSumm = 0
		#end
		SummFile = self.options['Summary file'][1]
		if (iSumm != 0):
			if os.path.isfile(SummFile):
				os.remove(SummFile)
			#end
			ierror = snopt.openunit(iSumm, numpy.array(SummFile), numpy.array('new'), numpy.array('sequential'))
			if (ierror != 0):
				raise IOError('Failed to properly open %s, ierror = %3d' %(SummFile,ierror))
			#end
		#end
		lencw = 500
		leniw = 500 + 100 * (ncon+nvar)
		lenrw = 500 + 200 * (ncon+nvar)
		
		self.options['Total integer workspace'][1] = leniw
		self.options['Total real workspace'][1] = lenrw
		
		string = []
		for i in range(lencw):
			string.append('        ')
		#end
		cw = numpy.array(string,'c')
		iw = numpy.zeros(leniw, 'i')
		rw = numpy.zeros(lenrw, numpy.float)
		snopt.sninit(iPrint, iSumm, cw, iw, rw)
		
		# Memory allocation
		ne = ncon*nvar
		nnCon = numpy.array(ncon, numpy.int)
		nnObj = numpy.array(nvar, numpy.int)
		nnJac = numpy.array(nvar, numpy.int)
		neGcon = nnCon*nnJac
		iExit = 0
		mincw, miniw, minrw,cw = snopt.snmemb(iExit, ncon, nvar, ne, neGcon, nnCon, nnJac, nnObj, cw, iw, rw)
		if ((minrw > lenrw) or (miniw > leniw) or (mincw > lencw)):
			if (mincw > lencw):
				lencw = mincw
				string = ""
				for i in range(lencw):
					string = string + "        "
				# end for
				cw = numpy.transpose(numpy.reshape(numpy.array(string),(lencw,8)))
			#end
			if (miniw > leniw):
				leniw = miniw
				iw = numpy.zeros(leniw, 'i')
			#end
			if (minrw > lenrw):
				lenrw = minrw        
				rw = numpy.zeros(lenrw, numpy.float)
			#end
			snopt.sninit(iPrint, iSumm, cw, iw, rw)
		#end
		
		
		# Set Options
		inform = numpy.array([-1], numpy.int)
		for i in xrange(len(self.set_options)):
			name = self.set_options[i][0]
			value = self.set_options[i][1]
			if isinstance(value, str):
				if (name == 'Start'):
					if (value == 'Cold'):
						snopt.snset('Cold start', iPrint, iSumm, inform, cw, iw, rw)
					elif (value == 'Warm'):
						snopt.snset('Warm start', iPrint, iSumm, inform, cw, iw, rw)
					#end
				elif (name == 'Problem Type'):
					if (value == 'Minimize'):
						snopt.snset('Minimize', iPrint, iSumm, inform, cw, iw, rw)
					elif (value == 'Maximize'):
						snopt.snset('Maximize', iPrint, iSumm, inform, cw, iw, rw)
					elif (value == 'Feasible point'):
						snopt.snset('Feasible point', iPrint, iSumm, inform, cw, iw, rw)
					#end
				elif (name == 'Print file'):
					snopt.snset(name+' '+'%d'%(iPrint), iPrint, iSumm, inform, cw, iw, rw)
				elif (name == 'Summary file'):
					snopt.snset(name+' '+'%d'%(iSumm), iPrint, iSumm, inform, cw, iw, rw)
				else:
					snopt.snset(name+' '+value, iPrint, iSumm, inform, cw, iw, rw)
				#end
			elif isinstance(value, float):
				snopt.snsetr(name, value, iPrint, iSumm, inform, cw, iw, rw)
			elif isinstance(value, int):
				if (name=='iPrint') and iPrint == 0:
					pass
				elif (name=='iSumm') and iSumm == 0:
					pass
				else:
					snopt.snseti(name, value, iPrint, iSumm, inform, cw, iw, rw)
				#end
			elif isinstance(value, type(None)):
				snopt.snset(name, iPrint, iSumm, inform, cw, iw, rw)
			#end
		#end
		
		
		# Setup argument list values
		start = numpy.array(self.options['Start'][1])
		nName = numpy.array([1], numpy.int)
		nnCon = numpy.array(ncon, numpy.int)
		nnObj = numpy.array(nvar, numpy.int)
		nnJac = numpy.array(nvar, numpy.int)
		iObj   = numpy.array([0], numpy.int)
		ObjAdd = numpy.array([0.], numpy.float)
		ProbNm = numpy.array(self.name)
		a = numpy.zeros(ne, numpy.float)
		ha = numpy.zeros(ne, 'i')
		ine = 0
		for j in range(nvar):
			for i in range(ncon):
				ha[ine] = i + 1
				ine = ine + 1
			#end            
		#end            
		ka = numpy.zeros(nvar+1, 'i')
		ka[0] = 1
		for i in range(1,nvar+1):
			ka[i] = ka[i-1] + ncon
		#end
		xs = numpy.concatenate((xs, numpy.zeros(ncon,numpy.float)))
		bl = numpy.concatenate((blx, blc))
		bu = numpy.concatenate((bux, buc))
		lencu = numpy.array([1], numpy.int)
		leniu = numpy.array([1], numpy.int)
		lenru = numpy.array([1], numpy.int)
		cu = numpy.array(["        "],'c')
		iu = numpy.zeros([leniu[0]], numpy.int)
		ru = numpy.zeros([lenru[0]], numpy.float)
		hs = numpy.zeros(nvar+ncon, 'i')
		Names = numpy.array(["        "],'c')
		pi = numpy.zeros(ncon, numpy.float)
		rc = numpy.zeros(nvar+ncon, numpy.float)
		#inform = numpy.array([-1], numpy.int)
		mincw = numpy.array([0], numpy.int)
		miniw = numpy.array([0], numpy.int)
		minrw = numpy.array([0], numpy.int)
		nS = numpy.array([0], numpy.int)
		ninf = numpy.array([0], numpy.int)
		sinf = numpy.array([0.], numpy.float)
		
		
		# Run SNOPT
		t0 = time.time()
		snopt.snoptc(start, nnCon, nnObj, nnJac, iObj, ObjAdd, ProbNm,
			snfuncgrag, a, ha, ka, bl, bu, Names, hs, xs, pi, rc, inform, 
			mincw, miniw, minrw, nS, ninf, sinf, ff, cu, iu, ru, cw, 
			iw, rw)
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
			
			if (iPrint != 0):
				snopt.closeunit(self.options['iPrint'][1])
			#end
			if (iSumm != 0):
				snopt.closeunit(self.options['iSumm'][1])
			#end
		#end
		
		# Store Results
		sol_inform = {}
		sol_inform['value'] = inform
		sol_inform['text'] = self.getInform(inform)
		
		if store_sol:
			
			sol_name = 'SNOPT Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_evals = 0
			
			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = xs[i]
				i += 1
			#end
			
			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = ff[i]
				i += 1
			#end
			
			sol_cons = copy.deepcopy(opt_problem._constraints)
			i = 0
			for key in sol_cons.keys():
				sol_cons[key].value = xs[nvar+i]
				i += 1
			#end
			
			if ncon > 0:
				sol_lambda = numpy.zeros(ncon ,float)
				for i in xrange(ncon):
					sol_lambda[i] = pi[i]
				#end
			else:
				sol_lambda = {}
			#end
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, Sensitivities=sens_type,
				myrank=myrank, arguments=args, **kwargs)
			
		#end
		
		return ff, xs[0:nvar], sol_inform
		
		
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 07, 2008 - Ruben E. Perez
		'''
		
		# 
		self.set_options.append([name,value])
		
		
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
		mjr_code = (infocode[0]/10)*10
		mnr_code = infocode[0] - 10*mjr_code
		try:
			inform_text = self.informs[mjr_code]
		except:
			inform_text = 'Unknown Exit Status'
		#end
		
		return inform_text
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		'''
		
		# 
		iPrint = self.options['iPrint'][1]
		iSumm = self.options['iSumm'][1]
		if (iPrint != 0):
			snopt.pyflush(iPrint)
		#end
		if (iSumm != 0):
			snopt.pyflush(iSumm)
		#end
		
		
	def SaveSpecsFile(self, filename):
		
		'''
		Save Options in a SPECS File
		
		Keyword arguments:
		-----------------
		filename -> STRING: Spec File Name
		
		Documentation last updated:  May. 08, 2008 - Ruben E. Perez
		'''
		
		# 
		
		
		
	def ReadSpecsFile(self, filename):
		
		'''
		Read in a SPECS file
		
		Keyword arguments:
		-----------------
		filename -> STRING: Spec File Name
		
		Documentation last updated:  May. 08, 2008 - Ruben E. Perez
		'''
		
		## 
		#if (os.path.isfile(filename)):
		#	
		#	try :
		#		
		#	except:
		#		raise IOError('Failed to properly open %s, ierror = %3d' %(filename,ierror))
		#	#end
		#	
		#	print 'Options had been updated based on %s file\n' %(filename)
		#	
		#else:
		#	raise IOError('Error: could not find file %s' %(filename))
		##end
	


#==============================================================================
# SNOPT Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test SNOPT
	print 'Testing ...'
	snopt = SNOPT()
	print snopt
	
