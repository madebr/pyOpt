#!/usr/bin/env python
'''
pyALPSO - A Python pyOpt interface to ALPSO.

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 16/05/2008 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter Jansen (PJ)

History
-------
	v. 1.0 	- Initial Class Creation (RP, 2008)
	v. 1.1	- Integrate to pyOpt Framework (RP, 2008)
'''

__version__ = '$Revision: $'

'''
To Do:
	- 
'''

# =============================================================================
# ALPSO Library
# =============================================================================
#try:
#	import alpso
#except:
#	pass
##end
#try:
#	import alpso_spm
#except:
#	pass
##end
#try:
#	import alpso_dpm
#except:
#	pass
##end

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
# ALPSO Optimizer Class
# =============================================================================
class ALPSO(Optimizer):
	
	'''
	ALPSO Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		ALPSO Optimizer Class Initialization
		
		**Keyword arguments:**
		
		- pll_type -> STR: ALPSO Parallel Implementation (None, SPM- Static, DPM- Dynamic, POA-Parallel Analysis), *Default* = None
		
		Documentation last updated:  February. 2, 2011 - Ruben E. Perez
		'''
		
		#
		if (pll_type == None):
			
			try:
				import alpso as alpso
			except:
				raise ImportError('pyALPSO: ALPSO shared library failed to import')
			#end
			
			name = 'ALPSO'
			self.alpso = alpso
		elif (pll_type.upper() == 'SPM'):
			
			try:
				import alpso_spm
				from mpi4py import MPI
			except:
				raise ImportError('pyALPSO: ALPSO SPM shared library failed to import')
			#end
			
			name = 'ALPSO - SPM'
			self.alpso = alpso_spm
		elif (pll_type.upper() == 'DPM'):
			
			#if not 'alpso_dpm' in sys.modules:
			#	raise ImportError('pyALPSO: ALPSO DPM shared library failed to import')
			##end
			
			try:
				import alpso_dpm
				from mpi4py import MPI
			except:
				raise ImportError('pyALPSO: ALPSO DPM shared library failed to import')
			#end
			
			name = 'ALPSO - DPM'
			self.alpso = alpso_dpm
		elif (pll_type.upper() == 'POA'):
			
			try:
				import alpso_poa
				from mpi4py import MPI
			except:
				raise ImportError('pyALPSO: ALPSO POA shared library failed to import')
			#end
			
			name = 'ALPSO - POA'
			self.alpso = alpso_poa
		else:
			raise ValueError("pll_type must be either None,'SPM', 'DPM' or 'POA'")
		#end
		
		category = 'Global Optimizer'
		def_opts = {
		'SwarmSize':[int,40],			# Number of Particles (Depends on Problem dimensions)
		'maxOuterIter':[int,200],		# Maximum Number of Outer Loop Iterations (Major Iterations)
		'maxInnerIter':[int,6],			# Maximum Number of Inner Loop Iterations (Minor Iterations)
		'minInnerIter':[int,6],			# Minimum Number of Inner Loop Iterations (Dynamic Inner Iterations)
		'dynInnerIter':[int,0],			# Dynamic Number of Inner Iterations Flag
		'stopCriteria':[int,1],			# Stopping Criteria Flag (0 - maxIters, 1 - convergence)
		'stopIters':[int,5],			# Consecutively Number of Iterations for which the Stopping Criteria must be Satisfied
		'etol':[float,1e-3],			# Absolute Tolerance for Equality constraints 
		'itol':[float,1e-3],			# Absolute Tolerance for Inequality constraints 
		#'ltol':[float,1e-2],			# Absolute Tolerance for Lagrange Multipliers
		'rtol':[float,1e-2],			# Relative Tolerance for Lagrange Multipliers
		'atol':[float,1e-2],			# Absolute Tolerance for Lagrange Function
		'dtol':[float,1e-1],			# Relative Tolerance in Distance of All Particles to Terminate (GCPSO)
		'printOuterIters':[int,0],		# Number of Iterations Before Print Outer Loop Information
		'printInnerIters':[int,0],		# Number of Iterations Before Print Inner Loop Information
		'rinit':[float,1.0],			# Initial Penalty Factor
		'xinit':[int,0],				# Initial Position Flag (0 - no position, 1 - position given)
		'vinit':[float,1.0],			# Initial Velocity of Particles in Normalized [-1,1] Design Space
		'vmax':[float,2.0],				# Maximum Velocity of Particles in Normalized [-1,1] Design Space
		'c1':[float,2.0],				# Cognitive Parameter
		'c2':[float,1.0],				# Social Parameter
		'w1':[float,0.99],				# Initial Inertia Weight
		'w2':[float,0.55],				# Final Inertia Weight
		'ns':[int,15],					# Number of Consecutive Successes in Finding New Best Position of Best Particle Before Search Radius will be Increased (GCPSO)
		'nf':[int,5],					# Number of Consecutive Failures in Finding New Best Position of Best Particle Before Search Radius will be Increased (GCPSO)
		'dt':[float,1.0],				# Time step
		'vcrazy':[float,1e-4],			# Craziness Velocity (Added to Particle Velocity After Updating the Penalty Factors and Langangian Multipliers)
		'fileout':[int,1],				# Flag to Turn On Output to filename
		'filename':[str,'ALPSO.out'],	# We could probably remove fileout flag if filename or fileinstance is given
		'seed':[float,0],				# Random Number Seed (0 - Auto-Seed based on time clock)
		'HoodSize':[int,40],			# Number of Neighbours of Each Particle
		'HoodModel':[str,'gbest'],		# Neighbourhood Model (dl/slring - Double/Single Link Ring, wheel - Wheel, Spatial - based on spatial distance, sfrac - Spatial Fraction)
		'HoodSelf':[int,1],				# Selfless Neighbourhood Model (0 - Include Particle i in NH i, 1 - Don't Include Particle i)
		'Scaling':[int,1],				# Design Variables Scaling Flag (0 - no scaling, 1 - scaling between [-1,1]) 
		}
		informs = {}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		#
		if (self.name in ('ALPSO - SPM','ALPSO - DPM','ALPSO - POA')):
			self.myrank = MPI.COMM_WORLD.Get_rank()
		else:
			self.myrank = 0
		#end
		
		
	def __solve__(self, opt_problem={}, store_sol=True, disp_opts=False, xstart=[], store_hst=False, hot_start=False, *args, **kwargs):
		
		'''
		Run Optimizer (Optimize Routine)
		
		**Keyword arguments:**
		
		- opt_problem -> INST: Optimization instance
		- store_sol -> BOOL: Store solution in Optimization class flag, *Default* = True 
		- disp_opts -> BOOL: Flag to display options in solution text, *Default* = False
		- xstart ->  :  , *Default* = []
		- store_hst -> BOOL/STR: Flag/filename to store optimization history, *Default* = False
		- hot_start -> BOOL/STR: Flag/filename to read optimization history, *Default* = False
		
		Additional arguments and keyword arguments are passed to the objective function call.
		
		Documentation last updated:  February. 2, 2011 - Ruben E. Perez
		'''
		
		#
		if kwargs.has_key('display_opts'):
			sol_dispOpt = kwargs['display_opts']
			del kwargs['display_opts']
		else:
			sol_dispOpt = False
		#end
		
		myrank = self.myrank
		
		#
		def_fname = self.options['filename'][1].split('.out')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		
		#======================================================================
		# ALPSO - Objective/Constraint Values Function
		#======================================================================
		def objconfunc(x):
			
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
			
			# Evaluate User Function
			[ff,gg,fail] = opt_problem.obj_fun(xn, *args, **kwargs)
			
			# 
			g = numpy.zeros(len(opt_problem._constraints.keys()),float)
			if (fail == 1):
				# Objective Assigment
				f = inf
				# Constraints Assigment
				for i in xrange(len(opt_problem._constraints.keys())):
					g[i] = inf
				#end
			else:
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
			#end
			
			return f,g
		
		
		
		# Variables Handling
		n = len(opt_problem._variables.keys())
		xl = numpy.zeros(n,float)
		xu = numpy.zeros(n,float)
		type = numpy.zeros(n,int)
		i = 0
		for key in opt_problem._variables.keys():
			xl[i] = opt_problem._variables[key].lower
			xu[i] = opt_problem._variables[key].upper
			if opt_problem._variables[key].type == 'c':
				type[i] = 0
			else:
				type[i] = 1
			#end
			i += 1
		#end
		
		# Variables Groups Handling 
		if opt_problem.use_groups:
			group_ids = {}
			k = 0
			for key in opt_problem._vargroups.keys():
				group_len = len(opt_problem._vargroups[key]['ids'])
				group_ids[opt_problem._vargroups[key]['name']] = [k,k+group_len]
				k += group_len
			#end
		#end
		
		# Constraints Handling
		m = len(opt_problem._constraints.keys())
		me = 0
		#i = 0
		if m > 0:
			for key in opt_problem._constraints.keys():
				if opt_problem._constraints[key].type == 'e':
					me += 1
				#end
			#end
			#i += 1
		#end
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		
		
		# Setup argument list values
		nob = self.options['SwarmSize'][1]
		nhn = self.options['HoodSize'][1]
		if (self.options['HoodModel'][1].lower() in ['gbest','dlring','slring','wheel','spatial','sfrac']):
			nhm = self.options['HoodModel'][1].lower()
		else:
			raise IOError('Incorrect Neighboorhood Model Setting')
		#end
		nhs = self.options['HoodSelf'][1]
		imax = self.options['maxOuterIter'][1]
		cmax = self.options['maxInnerIter'][1]
		cmin = self.options['minInnerIter'][1]
		dyniI = self.options['dynInnerIter'][1]
		if (dyniI == 0):
			cmin = cmax
		#end
		if (self.options['stopCriteria'][1]>=0 and self.options['stopCriteria'][1]<=1):
			stop = self.options['stopCriteria'][1]
		else:
			raise IOError('Incorrect Stopping Criteria Setting')
		#end
		nstop = self.options['stopIters'][1]
		etol = self.options['etol'][1]
		itol = self.options['itol'][1]
		#ltol = self.options['ltol'][1]
		rtol = self.options['rtol'][1]
		atol = self.options['atol'][1]
		dtol = self.options['dtol'][1]
		oout = self.options['printOuterIters'][1]
		iout = self.options['printInnerIters'][1]
		rinit = self.options['rinit'][1]
		xinit = self.options['xinit'][1]
		vinit = self.options['vinit'][1]
		vmax = self.options['vmax'][1]
		c1 = self.options['c1'][1]
		c2 = self.options['c2'][1]
		w1 = self.options['w1'][1]
		w2 = self.options['w2'][1]
		ns = self.options['ns'][1]
		nf = self.options['nf'][1]
		vcrazy = self.options['vcrazy'][1]
		if (self.options['fileout'][1]>=0 and self.options['fileout'][1]<=3):
			fileout = self.options['fileout'][1]
		else:
			raise IOError('Incorrect fileout Setting')
		#end
		filename = self.options['filename'][1]
		
		seed = self.options['seed'][1]
		if seed == 0:
			seed = time.time()
		#end
		if self.h_start:
			seed = hos_file.read(-1,ident=['seed'])[0]['seed'][0][0]
		#end
		scale = self.options['Scaling'][1]
		#dt = self.options['dt'][1]
		xs = []
		if (xinit == 1):
			xst = []
			for key in opt_problem._variables.keys():
				xst.append(opt_problem._variables[key].value)
			#end
			xst = numpy.array(xst)
			if (xstart == []):
				xs = xst
			elif (len(xstart) < nob):
				if isinstance(xstart,list):
					xstart.append(xst)
				elif isinstance(xstart,numpy.ndarray):
					xstart = numpy.concatenate((xstart,xst.reshape(1,n)))
				#end
				xs = xstart
			else:
				xs = xstart
			#end
		#end
		
		
		# Run ALPSO
		t0 = time.time()
		opt_x,opt_f,opt_g,opt_lambda,nfevals,rseed = self.alpso.alpso(n,m,me,
			type,xs,xl,xu,nob,nhn,nhm,imax,cmax,cmin,stop,nstop,etol,itol,
			rtol,atol,dtol,oout,iout,rinit,vinit,vmax,c1,c2,w1,w2,ns,nf,vcrazy,
			fileout,filename,log_file,hos_file,seed,scale,nhs,objconfunc)
		sol_time = time.time() - t0 
		sol_evals = nfevals
		
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
		if store_sol:
			
			sol_name = 'ALPSO Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_inform = {}
			#sol_inform['value'] = inform
			#sol_inform['text'] = self.getInform(inform)
			
			
			sol_vars = copy.deepcopy(opt_problem._variables)
			i = 0
			for key in sol_vars.keys():
				sol_vars[key].value = opt_x[i]
				i += 1
			#end
			
			sol_objs = copy.deepcopy(opt_problem._objectives)
			i = 0
			for key in sol_objs.keys():
				sol_objs[key].value = opt_f	# Note: takes only one!
				i += 1
			#end
			
			if m > 0:
				sol_cons = copy.deepcopy(opt_problem._constraints)
				i = 0
				for key in sol_cons.keys():
					sol_cons[key].value = opt_g[i]
					i += 1
				#end
			else:
				sol_cons = {}
			#end
			
			if m > 0:
				sol_lambda = numpy.zeros(m ,float)
				for i in xrange(m):
					sol_lambda[i] = opt_lambda[i]
				#end
			else:
				sol_lambda = {}
			#end
			
			
			opt_problem.addSol(self.__class__.__name__, sol_name, objfunc, sol_time, 
				sol_evals, sol_inform, sol_vars, sol_objs, sol_cons, sol_options, 
				display_opts=disp_opts, Lambda=sol_lambda, Seed=rseed, myrank=self.myrank,
				arguments=args, **kwargs)
		#end
		
		
		return opt_f, opt_x, {'opt_g':opt_g,'fevals':sol_evals,'time':sol_time}
		
		
	def _on_setOption(self, name, value):
		
		'''
		Set Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 16, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getOption(self, name):
		
		'''
		Get Optimizer Option Value (Optimizer Specific Routine)
		
		Documentation last updated:  May. 16, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_getInform(self, infocode):
		
		'''
		Get Optimizer Result Information (Optimizer Specific Routine)
		
		Keyword arguments:
		-----------------
		id -> STRING: Option Name
		
		Documentation last updated:  May. 16, 2008 - Ruben E. Perez
		'''
		
		pass
		
		
	def _on_flushFiles(self):
		
		'''
		Flush the Output Files (Optimizer Specific Routine)
		
		Documentation last updated:  August. 09, 2009 - Ruben E. Perez
		'''
		
		pass
	


#==============================================================================
# ALPSO Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test ALPSO
	print 'Testing ...'
	alpso = ALPSO()
	print alpso
	
