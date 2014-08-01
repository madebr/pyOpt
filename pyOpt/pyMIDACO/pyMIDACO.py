#!/usr/bin/env python
'''
pyMIDACO - A Python pyOpt interface to MIDACO. 

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.3   $Date: 31/07/2014 21:00$


Tested on:
---------
Linux with gfortran
Linux with pathf95
Win32 with gfortran
Mac with g95

Developers:
-----------
- Dr. Ruben E. Perez (RP)

History
-------
	v. 1.0	- Initial Class Creation (RP, 2009)
	v. 1.1  - Updated Functionality for MIDACO v.3.0 (RP, 2012)
	v. 1.2  - Updated Functionality for MIDACO v.4.0 (RP, 2014)
	v. 1.3	- Unconstrained Problems Support (RP, 2014)
'''

__version__ = '$Revision: $'

'''
To Do:
	- implement function parallelization, i.e. L > 1
'''

# =============================================================================
# MIDACO Library
# =============================================================================
try:
	import midaco
except:
	raise ImportError('MIDACO shared library failed to import')
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
# MIDACO Optimizer Class
# =============================================================================
class MIDACO(Optimizer):
	
	'''
	MIDACO Optimizer Class - Inherited from Optimizer Abstract Class
	'''
	
	def __init__(self, pll_type=None, *args, **kwargs):
		
		'''
		MIDACO Optimizer Class Initialization
		
		**Keyword arguments:**
		
		- pll_type -> STR: Parallel Implementation (None, 'POA'-Parallel Objective Analysis), *Default* = None
		
		Documentation last updated:  Feb. 16, 2010 - Peter W. Jansen
		'''
		
		#
		if (pll_type == None):
			self.poa = False
			self.spm = False
		elif (pll_type.upper() == 'POA'):
			self.poa = True
			self.spm = False
		elif (pll_type.upper() == 'SPM'):
			self.poa = False
			self.spm = True
		else:
			raise ValueError("pll_type must be either None, 'POA' or 'SPM'")
		#end
		
		#
		name = 'MIDACO'
		category = 'Global Optimizer'
		def_opts = {
		# MIDACO Options
		'ACC':[float,0],       			    # Accuracy for constraint violation (0 - default)
		'ISEED':[int,0],            		# Seed for random number generator  (e.g. ISEED = 0,1,2,3,...)
		'FSTOP':[int,0],					# Objective Function Stopping Value (0 - disabled)
		'AUTOSTOP':[int,0],					# Automatic stopping criteria (0 - disable, 1 to 500 - from local to global)
		'ORACLE':[float,0],					# Oracle parameter for constrained problems (0 - Use internal default)
		'FOCUS':[int,0],					# Focus of MIDACO search process around best solution (0 - Use internal default)
		'ANTS':[int,0],						# Number of iterates (ants) per generation (0 - Use internal default)
		'KERNEL':[int,0],					# Size of the solution archive (0 - Use internal default)
		'CHARACTER':[int,0],				# Internal custom parameters (0 - Use internal default, 1 - IP problems, 2 - NLP problems, 3 - MINLP problems)
		'MAXEVAL':[int,10000],				# Maximum function evaluations
		'MAXTIME':[int,86400],				# Maximum time limit, in seconds
		'IPRINT':[int,1],					# Output Level (<0 - None, 0 - Screen, 1 - File(s))
		'PRINTEVAL':[int,1000],				# Print-Frequency for current best solution
		'IOUT1':[int,36],					# History output unit number
		'IOUT2':[int,37],					# Best solution output unit number
		'IFILE1':[str,'MIDACO_HIST.out'],	# History output file name
		'IFILE2':[str,'MIDACO_BEST.out'],	# Best output file name
		'LKEY':[str,'MIDACO_LIMITED_VERSION___[CREATIVE_COMMONS_BY-NC-ND_LICENSE]'],
		}
		informs = {
		1 : 'Feasible solution,   MIDACO was stopped by the user submitting ISTOP=1',
		2 : 'Infeasible solution, MIDACO was stopped by the user submitting ISTOP=1',
		3 : 'Feasible solution,   MIDACO stopped automatically using AUTOSTOP option',
		4 : 'Infeasible solution,   MIDACO stopped automatically using AUTOSTOP option',
		5 : 'Feasible solution,   MIDACO stopped automatically by FSTOP',
		51 : 'WARNING: Some X(i)  is greater/lower than +/- 1.0D+12 (try to avoid huge values!)',
		52 : 'WARNING: Some XL(i) is greater/lower than +/- 1.0D+12 (try to avoid huge values!)',
		53 : 'WARNING: Some XU(i) is greater/lower than +/- 1.0D+12 (try to avoid huge values!)',
		61 : 'WARNING: Some X(i)  should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)',
		62 : 'WARNING: Some XL(i) should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)',
		63 : 'WARNING: Some XU(i) should be discrete (e.g. 1.000) , but is continuous (e.g. 1.234)',
		71 : 'WARNING: Some XL(i) = XU(I) (fixed variable)',
		81 : 'WARNING: F(X) has value NaN for starting point X (sure your problem is correct?)',
		82 : 'WARNING: Some G(X) has value NaN for starting point X (sure your problem is correct?)',
		91 : 'WARNING: FSTOP is greater/lower than +/- 1.0D+8',
		92 : 'WARNING: ORACLE is greater/lower than +/- 1.0D+8',
		101 : 'ERROR: L    <= 0 or L > 1.0D+6',
		102 : 'ERROR: N    <= 0 or N > 1.0D+6',
		103 : 'ERROR: NINT <  0',
		104 : 'ERROR: NINT >  N',
		105 : 'ERROR: M    <  0 or M > 1.0D+6',
		106 : 'ERROR: ME   <  0',
		107 : 'ERROR: ME   >  M',
		201 : 'ERROR: some X(i)  has type NaN',
		202 : 'ERROR: some XL(i) has type NaN',
		203 : 'ERROR: some XU(i) has type NaN',
		204 : 'ERROR: some X(i) < XL(i)',
		205 : 'ERROR: some X(i) > XU(i)',
		206 : 'ERROR: some XL(i) > XU(i)',
		301 : 'ERROR: ACC < 0   or   ACC > 1.0D+6',
		302 : 'ERROR: ISEED < 0   or   ISEED > 1.0D+12',
		303 : 'ERROR: FSTOP greater/lower than +/- 1.0D+12',
		304 : 'ERROR: AUTOSTOP < 0   or   AUTOSTOP > 1.0D+6',
		305 : 'ERROR: ORACLE greater/lower than +/- 1.0D+12',
		306 : 'ERROR: |FOCUS| < 1   or   FOCUS > 1.0D+12',
		307 : 'ERROR: ANTS < 0   or   ANTS > 1.0D+8',
		308 : 'ERROR: KERNEL < 0   or   KERNEL > 100',
		309 : 'ERROR: ANTS < KERNEL',
		310 : 'ERROR: ANTS > 0 but KERNEL = 0',
		311 : 'ERROR: KERNEL > 0 but ANTS = 0',
		312 : 'ERROR: CHARACTER < 0   or   CHARACTER > 1000',
		313 : 'ERROR: some MIDACO parameters has type NaN',
		401 : 'ERROR: ISTOP < 0 or ISTOP > 1',
		501 : 'ERROR: Double precision work space size LRW is too small (see below LRW), RW must be at least of size LRW = 200*N+2*M+1000',
		601 : 'ERROR: Integer work space size LIW is too small (see below LIW), IW must be at least of size LIW = 2*N+L+100',
		701 : 'ERROR: Input check failed! MIDACO must be called initially with IFAIL = 0',
		801 : 'ERROR: L > LMAX (user must specifiy LMAX below in the MIDACO source code)',
		802 : 'ERROR: L*M+1 > LXM (user must specifiy LXM below in the MIDACO source code)',
		900 : 'ERROR: Invalid or corrupted LICENSE_KEY',
		999 : 'ERROR: N > 4. The free test version is limited up to 4 variables',
		}
		Optimizer.__init__(self, name, category, def_opts, informs, *args, **kwargs)
		
		
	def __solve__(self, opt_problem={}, store_sol=True, disp_opts=False, store_hst=False, hot_start=False, *args, **kwargs):
		
		'''
		Run Optimizer (Optimize Routine)
		
		**Keyword arguments:**
		
		- opt_problem -> INST: Optimization instance
		- store_sol -> BOOL: Store solution in Optimization class flag, *Default* = True 
		- disp_opts -> BOOL: Flag to display options in solution text, *Default* = False
		- store_hst -> BOOL/STR: Flag/filename to store optimization history, *Default* = False
		- hot_start -> BOOL/STR: Flag/filename to read optimization history, *Default* = False
		
		Additional arguments and keyword arguments are passed to the objective function call.
		
		Documentation last updated:  February. 17, 2011 - Peter W. Jansen
		'''
		#
		if self.poa or self.spm:
			try:
				import mpi4py
				from mpi4py import MPI
			except ImportError:
				print 'pyMIDACO: Parallel objective Function Analysis requires mpi4py'
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
			self.myrank = comm.Get_rank()
			if self.poa:
				nproc = 1
			else:
				nproc = comm.Get_size()
			#end
		else:
			self.pll = False
			nproc = 1
			self.myrank = 0
		#end
		
		myrank = self.myrank
		
		# 
		def_fname = self.options['IFILE1'][1].split('.')[0]
		hos_file, log_file, tmp_file = self._setHistory(opt_problem.name, store_hst, hot_start, def_fname)
		
		
		#======================================================================
		# MIDACO - Objective/Constraint Values Function
		#======================================================================
		def objfun(l,n,m,x,f,g):
			
			x = numpy.reshape(x,(l,-1))
			f = numpy.reshape(f,(l,-1))
			g = numpy.reshape(g,(l,-1))
			
			# Variables Groups Handling
			if not self.poa:
				mxi = myrank
			else:
				mxi = 0
			#end
			if opt_problem.use_groups:
				xg = {}
				for group in group_ids.keys():
					if (group_ids[group][1]-group_ids[group][0] == 1):
						xg[group] = x[mxi,group_ids[group][0]]
					else:
						xg[group] = x[mxi,group_ids[group][0]:group_ids[group][1]]
					#end
				#end
				xn = xg
			else:
				xn = x[mxi,:]
			#end
			
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
							g[proc,:] = gg[:]
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
				
				# 
				if (fail == 1):
					# Objective Assigment
					f[mxi] = inf
					# Constraints Assigment (negative gg as midaco uses g(x) >= 0)
					for i in xrange(len(opt_problem._constraints.keys())):
						g[mxi,i] = -inf
					#end
				else:
					# Objective Assigment
					if isinstance(ff,complex):
						f[mxi] = ff.astype(float)
					else:
						f[mxi] = ff
					#end
					# Constraints Assigment (negative gg as midaco uses g(x) >= 0)
					for i in xrange(len(opt_problem._constraints.keys())):
						if isinstance(gg[i],complex):
							g[mxi,i] = -gg[i].astype(float)
						else:
							g[mxi,i] = -gg[i]
						#end
					#end
				#end
				
				if self.spm:
					send_buf = {}
					send_buf[myrank] = {'fi':f[mxi],'gi':g[mxi]}
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
								g[i,:] = p_results[proc][i]['gi']
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
						log_file.write(x[proc],'x')
						log_file.write(f[proc],'obj')
						log_file.write(g[proc],'con')
						log_file.write(fail,'fail')
					#end
				#end
			#end
			
			f = numpy.reshape(f,l)
			g = numpy.reshape(g,l*m)
			
			return f,g
		
		
		
		# Variables Handling
		nvar = len(opt_problem._variables.keys())
		nint = 0
		xl = []
		xu = []
		xx = []
		for key in opt_problem._variables.keys():
			if (opt_problem._variables[key].type == 'i'):
				nint += 1
			#end
			xl.append(opt_problem._variables[key].lower)
			xu.append(opt_problem._variables[key].upper)
			xx.append(opt_problem._variables[key].value)
		#end
		xl = numpy.array(xl)
		xu = numpy.array(xu)
		xx = numpy.array(xx*nproc)
		
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
		else:
			ncon = 1
			gg.append(0.0)
		#end
		gg = numpy.array(gg*nproc)
		
		# Objective Handling
		objfunc = opt_problem.obj_fun
		nobj = len(opt_problem._objectives.keys())
		
		ff = []
		for key in opt_problem._objectives.keys():
			ff.append(opt_problem._objectives[key].value)
		#end
		ff = numpy.array(ff*nproc)
		
		
		# Setup argument list values
		ll = numpy.array([nproc], numpy.int)
		nn = numpy.array([nvar], numpy.int)
		ni = numpy.array([nint], numpy.int)
		mm = numpy.array([ncon], numpy.int)
		me = numpy.array([neqc], numpy.int)
		opts = [0]*9
		opts[0] = self.options['ACC'][1]
		opts[1] = self.options['ISEED'][1]
		opts[2] = self.options['FSTOP'][1]
		opts[3] = self.options['AUTOSTOP'][1]
		opts[4] = self.options['ORACLE'][1]
		opts[5] = self.options['FOCUS'][1]
		opts[6] = self.options['ANTS'][1]
		opts[7] = self.options['KERNEL'][1]
		opts[8] = self.options['CHARACTER'][1]
		param = numpy.array([opts], numpy.float)
		maxeval = numpy.array([self.options['MAXEVAL'][1]], numpy.int)
		maxtime = numpy.array([self.options['MAXTIME'][1]], numpy.int)
		if (myrank == 0):
			iprint = numpy.array([self.options['IPRINT'][1]], numpy.int)
		else:
			iprint = numpy.array([-1], numpy.int)
		#end
		ifail = numpy.array([0], numpy.int)
		neval = numpy.array([0], numpy.int)
		if (self.options['PRINTEVAL'][1] > 0):
			printeval = numpy.array([self.options['PRINTEVAL'][1]], numpy.int)
		else:
			raise IOError('Incorrect PRINTEVAL Setting')
		#end
		iout1 = numpy.array([self.options['IOUT1'][1]], numpy.int)
		iout2 = numpy.array([self.options['IOUT2'][1]], numpy.int)
		ifile1 = self.options['IFILE1'][1]
		ifile2 = self.options['IFILE2'][1]
		if (iprint > 0):
			if os.path.isfile(ifile1):
				os.remove(ifile1)
			#end
			if os.path.isfile(ifile2):
				os.remove(ifile2)
			#end
		#end
		lkey = self.options['LKEY'][1]
		liw0 = 2*nn + ll + 1000
		liw = numpy.array([liw0], numpy.int)
		iw = numpy.zeros([liw], numpy.int)
		lrw0 = 200*nn + 2*mm + 1000
		lrw = numpy.array([lrw0], numpy.int)
		rw = numpy.zeros([lrw], numpy.float)
		
		
		# Run MIDACO
		t0 = time.time()
		midaco.midaco_wrap(ll,nn,ni,mm,me,xx,xl,xu,ff,gg,param,maxeval,maxtime,ifail,neval,iprint,printeval,iout1,iout2,ifile1,ifile2,lkey,liw,iw,lrw,rw,objfun)
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
			midaco.closeunit(self.options['IOUT1'][1])
			midaco.closeunit(self.options['IOUT2'][1])
		#end
		
		
		# Store Results
		if store_sol:
			
			sol_name = 'MIDACO Solution to ' + opt_problem.name
			
			sol_options = copy.copy(self.options)
			if sol_options.has_key('defaults'):
				del sol_options['defaults']
			#end
			
			sol_inform = {}
			sol_inform['value'] = ifail[0]
			sol_inform['text'] = self.getInform(ifail[0])
			
			sol_evals = neval
			
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
				display_opts=disp_opts, Lambda=sol_lambda, 
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
		iprint = self.options['IPRINT'][1]
		if (iprint > 0):
			midaco.pyflush(self.options['IOUT1'][1])
			midaco.pyflush(self.options['IOUT2'][1])
		#end
	


#==============================================================================
# MIDACO Optimizer Test
#==============================================================================
if __name__ == '__main__':
	
	# Test MIDACO
	print 'Testing ...'
	midaco = MIDACO()
	print midaco
	
