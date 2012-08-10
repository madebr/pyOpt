#!/usr/bin/env python
'''
alhso - Python Version of the Augmented Lagrangian Harmony Search Optimizer

hso if a global optimizer which solves problems of the form:

			min F(x)
				
	subject to: Gi(x)  = 0, i = 1(1)ME
				Gj(x) <= 0, j = ME+1(1)M
				xLB <= x <= xUB

Copyright (c) 2009-2011 by Dr. Ruben E. Perez
All rights reserved.
Revision: 1.2   $Date: 17/02/2009 21:00$


Tested on:
---------
- 

Developers:
-----------
- Dr. Ruben E. Perez (RP)
- 

History
-------
	v. 1.0	- Initial Code Development (RP, 2008)
	v. 1.1  - Extend Code to Handle Constraints (RP, 2008)
	v. 1.2  - Implemented Augmented Lagrangian Code (RP,2009)
		- Removed lambda convergance citeria (RP, 2009)
'''

__version__ = '$Revision: $'

'''
To Do:
	- Fix integer limits bug
	- How to handle tolerances definition change when scaled?
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys, random, time
import pdb
from math import floor

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================


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


#==============================================================================
# alhso function
#==============================================================================
def alhso(dimensions,constraints,neqcons,xtype,x0,xmin,xmax,
	memsize,maxoutiter,maxinniter,stopcriteria,stopiters,etol,
	itol,atol,rtol,prtoutiter,prtinniter,r0,hmcr,par,bw,
	fileout,filename,rseed,scale,objfunc):
	
	'''
	Python Version of the Augmented Lagrangian Harmony Search Optimizer
	
	Documentation last updated:  February. 19, 2009 - Ruben E. Perez
	'''
	
	# Set random number seed
	rand = random.Random()
	if rseed == {}:	
		rseed = time.time()
	#end
	
	rand.seed(rseed)
	
	# 
	if (fileout == 1):
		if (filename == ''):
			filename = 'Print.out'
		#end
		ofile = open(filename,'w')
	#end
	
	# 
	if (scale == 1):
		dbw = (xmax - xmin)/bw
		space_centre = numpy.zeros(dimensions,float)
		space_halflen = numpy.zeros(dimensions,float)
		for j in xrange(dimensions):
			space_centre[j] = (xmin[j] + xmax[j])/2.0
			space_halflen[j] = ((xmax[j] - xmin[j])/2.0)
		#end
		xmin = -numpy.ones(dimensions,float)
		xmax =  numpy.ones(dimensions,float)
		bw = (xmax - xmin)/dbw
	#end
	
	# Initialize Augmented Lagrange
	rp_val = numpy.ones(constraints, float)*r0
	lambda_val = numpy.zeros(constraints, float)
	lambda_old = numpy.zeros(constraints, float)
	
	# Initialize Harmony Memory
	HM = numpy.zeros((memsize,dimensions+1), float)
	discrete_i = []
	for i in xrange(memsize):
		for j in xrange(dimensions):
			HM[i,j] = xmin[j] + rand.random()*(xmax[j]-xmin[j])
			if (xtype[j] == 1):
				discrete_i.append(j)
			#end
		#end
	#end
	if (x0 != []):
		if (scale == 1):
			HM[:,:-1] = (x0[:] - space_centre)/space_halflen
		else:
			HM[:,:-1] = x0[:]
		#end
	#end
	
	# Initialize Harmony Memory Augmented Lagrange
	x_val = numpy.zeros(dimensions, float)
	x_tmp = numpy.zeros(dimensions, float)
	tau_val = numpy.zeros(constraints, float)
	nfevals = 0
	#best_L_val = 0
	for i in xrange(memsize):
		
		# Evaluate Ojective Function
		if (scale == 1):
			x_tmp = (HM[i,:-1] * space_halflen) + space_centre
		else:
			x_tmp = HM[i,:-1]
		#end
		for m in discrete_i:
			x_tmp[m] = floor(x_tmp[m] + 0.5)
		#end
		[f_val,g_val] = objfunc(x_tmp)
		nfevals += 1
		
		# Augmented Lagrangian Value
		L_val = f_val
		if (constraints > 0):
			
			# Equality Constraints
			for l in xrange(neqcons):
				tau_val[l] = g_val[l]
			#end
			
			# Inequality Constraints
			for l in xrange(neqcons,constraints):
				if (rp_val[l] != 0):
					if (g_val[l] > -lambda_val[l]/(2*rp_val[l])):
						tau_val[l] = g_val[l]
					else:
						tau_val[l] = -lambda_val[l]/(2*rp_val[l])
					#end
				else:
					tau_val[l] = g_val[l]
				#end
			#end
			
			#
			for l in xrange(constraints):
				L_val += lambda_val[l]*tau_val[l] + rp_val[l]*tau_val[l]**2
			#end
			
		#end
		
		# 
		HM[i,dimensions] = L_val
		
	#end
	
	# Initialize Best
	best_x_val = numpy.zeros(dimensions, float)
	best_f_val = []
	best_g_val = numpy.zeros(constraints, float)
	
	best_x_old = numpy.zeros(dimensions, float)
	best_f_old = []
	best_g_old = numpy.zeros(constraints, float)
	
	
	# Outer Optimization Loop
	k_out = 0
	kobj = 0
	iobj = 0
	stop_main_flag = 0
	while ((k_out < maxoutiter) and (stop_main_flag == 0)):
		
		k_out += 1
		
		# Inner Optimization Loop
		k_inn = 0
		while (k_inn < maxinniter):
			
			k_inn += 1
			
			# New Harmony Improvisation
			for j in xrange(dimensions):
				
				if ((rand.random() < hmcr) or (x0 != [] and k_out == 1)):
					
					# Harmony Memory Considering
					x_val[j] = HM[int(memsize*rand.random()),j]
					
					# Pitch Adjusting
					if (rand.random() <= par):
						if (rand.random() > 0.5):
							x_val[j] += rand.random()*bw[j]
						else:
							x_val[j] -= rand.random()*bw[j]
						#end
					#end
					
				else:
					
					# Random Searching
					x_val[j] = xmin[j] + rand.random()*(xmax[j]-xmin[j])
					
				#end
				
				# Check for improvisations out of range
				if (x_val[j] > xmax[j]): 
					x_val[j] = xmax[j]
				elif (x_val[j] < xmin[j]):
					x_val[j] = xmin[j]
				#end
				
			#end
			
			# Evaluate 
			if (scale == 1):
				x_tmp = (x_val * space_halflen) + space_centre
			else:
				x_tmp = x_val
			#end
			for m in discrete_i:
				x_tmp[m] = floor(x_tmp[m] + 0.5)
			#end
			[f_val,g_val] = objfunc(x_tmp)
			nfevals += 1
			
			# Lagrangian Value
			L_val = f_val
			if (constraints > 0):
				
				# Equality Constraints
				for l in xrange(neqcons):
					tau_val[l] = g_val[l]
				#end
				
				# Inequality Constraints
				for l in xrange(neqcons,constraints):
					if (rp_val[l] != 0):
						if (g_val[l] > -lambda_val[l]/(2*rp_val[l])):
							tau_val[l] = g_val[l]
						else:
							tau_val[l] = -lambda_val[l]/(2*rp_val[l])
						#end
					else:
						tau_val[l] = g_val[l]
					#end
				#end
				
				# 
				for l in xrange(constraints):
					L_val += lambda_val[l]*tau_val[l] + rp_val[l]*tau_val[l]**2
				#end
				
			#end
			
			# 
			feasible = True
			if (constraints > 0):
				for l in xrange(constraints):
					if (l < neqcons):
						if (abs(g_val[l]) > etol):
							feasible = False
							break
						#end
					else:
						if (g_val[l] > itol):
							feasible = False
							break
						#end
					#end
				#end
			#end
			
			# 
			if (feasible or (k_out == 1 and x0 != [])):
				
				# Harmony Memory Update
				hmax_num = 0
				hmax = HM[0,dimensions]
				for i in xrange(memsize):
					if (HM[i,dimensions] > hmax):
						hmax_num = i
						hmax = HM[i,dimensions]
					#end
				#end
				
				if (L_val < hmax):
					for j in xrange(dimensions):
						HM[hmax_num,j] = x_val[j]
					#end
					HM[hmax_num,dimensions] = L_val
				#end
				
				hmin_num = 0
				hmin = HM[0,dimensions]
				for i in xrange(memsize):
					if (HM[i,dimensions] < hmin):
						hmin_num = i
						hmin = HM[i,dimensions]
					#end
				#end
				
				if (L_val == hmin):
					
					best_x_val = x_val
					best_f_val = f_val
					best_g_val = g_val
					
					# Print Inner
					if (prtinniter != 0):
						# output to screen
						print '%d Inner Iteration of %d Outer Iteration' %(k_inn,k_out)
						print L_val
						
						if (scale == 1):
							x_tmp = (x_val * space_halflen) + space_centre
						else:
							x_tmp = x_val
						#end
						for m in discrete_i:
							x_tmp[m] = floor(x_tmp[m] + 0.5)
						#end
						print x_tmp
						
						print f_val
						print g_val
						print nfevals
					#end
					if (fileout == 1):
						# output to filename
						pass
					#end
					
					break
					
				#end
				
			#end
			
		#end
		
		# 
		if (best_f_val == [] and k_out == 1 and x0 == []):
			
			# Re-Initialize Harmony Memory
			HM = numpy.zeros((memsize,dimensions+1), float)
			for i in xrange(memsize):
				for j in xrange(dimensions):
					HM[i,j] = xmin[j] + rand.random()*(xmax[j]-xmin[j])
				#end
			#end
			
			# Re-Initialize Harmony Memory Augmented Lagrange
			for i in xrange(memsize):
				
				# Evaluate Ojective Function
				if (scale == 1):
					x_tmp = (HM[i,:-1] * space_halflen) + space_centre
				else:
					x_tmp = HM[i,:-1]
				#end
				for m in discrete_i:
					x_tmp[m] = floor(x_tmp[m] + 0.5)
				#end
				[f_val,g_val] = objfunc(x_tmp)
				nfevals += 1
				
				# Augmented Lagrangian Value
				L_val = f_val
				if (constraints > 0):
					
					# Equality Constraints
					for l in xrange(neqcons):
						tau_val[l] = g_val[l]
					#end
					
					# Inequality Constraints
					for l in xrange(neqcons,constraints):
						if (rp_val[l] != 0):
							if (g_val[l] > -lambda_val[l]/(2*rp_val[l])):
								tau_val[l] = g_val[l]
							else:
								tau_val[l] = -lambda_val[l]/(2*rp_val[l])
							#end
						else:
							tau_val[l] = g_val[l]
						#end
					#end
					
					#
					for l in xrange(constraints):
						L_val += lambda_val[l]*tau_val[l] + rp_val[l]*tau_val[l]**2
					#end
					
				#end
				
				# 
				HM[i,dimensions] = L_val
				
			#end
			
			# 
			k_out -= 1
			continue
			
		#end
		
		
		# Print Outer
		if (prtoutiter != 0 and numpy.mod(k_out,prtoutiter) == 0):
			
			# Output to screen
			print("="*80 + "\n")
			print("NUMBER OF ITERATIONS: %d\n" %(k_out))
			print("NUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" %(nfevals))
			print("OBJECTIVE FUNCTION VALUE:")
			print("\tF = %g\n" %(best_f_val))
			if (constraints > 0):
				# Equality Constraints
				print("EQUALITY CONSTRAINTS VALUES:")
				for l in xrange(neqcons):
					print("\tG(%d) = %g" %(l,best_g_val[l]))
				#end
				# Inequality Constraints
				print("\nINEQUALITY CONSTRAINTS VALUES:")
				for l in xrange(neqcons,constraints):
					print("\tH(%d) = %g" %(l,best_g_val[l]))
				#end
			#end
			print("\nLAGRANGIAN MULTIPLIERS VALUES:")
			for l in xrange(constraints):
				print("\tL(%d) = %g" %(l,lambda_val[l]))
			#end
			
			print("\nDESIGN VARIABLES VALUES:")
			if (scale == 1):
				x_tmp = (best_x_val[:] * space_halflen) + space_centre
			else:
				x_tmp = best_x_val[:]
			#end
			for m in discrete_i:
				x_tmp[m] = floor(x_tmp[m]+0.5)
			#end
			text = ''
			for j in xrange(dimensions):
				text += ("\tP(%d) = %9.3e\t" %(j,x_tmp[j]))
				if (numpy.mod(j+1,3) == 0):
					text +=("\n")
				#end
			#end
			print text
			print("="*80 + "\n")
		#end
		if (fileout == 1):
			# Output to filename
			ofile.write("\n" + "="*80 + "\n")
			ofile.write("\nNUMBER OF ITERATIONS: %d\n" %(k_out))
			ofile.write("\nNUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" %(nfevals))
			ofile.write("\nOBJECTIVE FUNCTION VALUE:\n")
			ofile.write("\tF = %g\n" %(best_f_val))
			if (constraints > 0):
				# Equality Constraints
				ofile.write("\nEQUALITY CONSTRAINTS VALUES:\n")
				for l in xrange(neqcons):
					ofile.write("\tG(%d) = %g\n" %(l,best_g_val[l]))
				#end
				# Inequality Constraints
				ofile.write("\nINEQUALITY CONSTRAINTS VALUES:\n")
				for l in xrange(neqcons,constraints):
					ofile.write("\tH(%d) = %g\n" %(l,best_g_val[l]))
				#end
			#end
			ofile.write("\nLAGRANGIAN MULTIPLIERS VALUES:\n")
			for l in xrange(constraints):
				ofile.write("\tL(%d) = %g\n" %(l,lambda_val[l]))
			#end
			
			ofile.write("\nDESIGN VARIABLES VALUES:\n")
			if (scale == 1):
				x_tmp = (best_x_val[:] * space_halflen) + space_centre
			else:
				x_tmp = best_x_val[:]
			#end
			for m in discrete_i:
				x_tmp[m] = floor(x_tmp[m]+0.5)
			#end
			text = ''
			for j in xrange(dimensions):
				text += ("\tP(%d) = %9.3e\t" %(j,x_tmp[j]))
				if (numpy.mod(j+1,3) == 0):
					text +=("\n")
				#end
			#end
			ofile.write(text)
			ofile.write("\n" + "="*80 + "\n")
			ofile.flush()
		#end
		
		
		# Test Constraint convergence
		stop_constraints_flag = 0
		if (constraints == 0):
			stop_constraints_flag = 1
		else:
			for l in xrange(neqcons):
				if (abs(best_g_val[l]) <= etol):
					stop_constraints_flag += 1 
				#end
			#end
			for l in xrange(neqcons,constraints):
				if (best_g_val[l] <= itol):
					stop_constraints_flag += 1
				#end
			#end
			if (stop_constraints_flag == constraints):
				stop_constraints_flag = 1
			else:
				stop_constraints_flag = 0
			#end
		#end
		
		# Test Position and Function convergence
		if (best_f_old == []):
			best_f_old = best_f_val
		#end
		stop_criteria_flag = 0
		if (stopcriteria == 1):
			
			# Absolute Change in Objective
			absfdiff = abs(best_f_val - best_f_old)
			if (absfdiff <= atol):
				kobj += 1
			else:
				kobj = 0
			#end
			
			# Relative Change in Objective
			if (abs(best_f_old) > 1e-10):
				if (abs(absfdiff/abs(best_f_old)) <= rtol):
					iobj += 1
				else:
					iobj = 0
				#end
			#end
			
			# 
			best_f_old = best_f_val
			
			# 
			if (kobj > stopiters or iobj > stopiters):
				stop_criteria_flag = 1
			else:
				stop_criteria_flag = 0
			#end
			
		#end
		
		# Test Convergence
		if (stop_constraints_flag == 1 and stop_criteria_flag == 1):
			stop_main_flag = 1
		else:
			stop_main_flag = 0
		#end
		
		
		# Update Augmented Lagrangian Terms
		if (stop_main_flag == 0):
			
			if (constraints > 0):
				
				# Tau for Best 
				for l in xrange(neqcons):
					tau_val[l] = best_g_val[l]
				#end
				for l in xrange(neqcons,constraints):
					if (best_g_val[l] > -lambda_val[l]/(2*rp_val[l])):
						tau_val[l] = best_g_val[l]
					else:
						tau_val[l] = -lambda_val[l]/(2*rp_val[l])
					#end
				#end
				
				# Update Lagrange Multiplier
				for l in xrange(constraints):
					lambda_old[l] = lambda_val[l]
					lambda_val[l] += 2*rp_val[l]*tau_val[l]
				#end
				
				# Update Penalty Factor
				for l in xrange(neqcons):
					if (abs(best_g_val[l]) > abs(best_g_old[l]) and abs(best_g_val[l]) > etol):
						rp_val[l] = 2.0*rp_val[l]
					elif (abs(best_g_val[l]) <= etol):
						rp_val[l] = 0.5*rp_val[l]
					#end
				#end
				for l in xrange(neqcons,constraints):
					if (best_g_val[l] > best_g_old[l] and best_g_val[l] > itol):
						rp_val[l] = 2.0*rp_val[l]
					elif (best_g_val[l] <= itol):
						rp_val[l] = 0.5*rp_val[l]
					#end
				#end
				
				# Apply Lower Bounds on rp
				for l in xrange(neqcons):
					if (rp_val[l] < 0.5*(abs(lambda_val[l])/etol)**0.5):
						rp_val[l] = 0.5*(abs(lambda_val[l])/etol)**0.5
					#end
				#end	
				for l in xrange(neqcons,constraints):
					if (rp_val[l] < 0.5*(abs(lambda_val[l])/itol)**0.5):
						rp_val[l] = 0.5*(abs(lambda_val[l])/itol)**0.5
					#end
				#end
				for l in xrange(constraints):
					if (rp_val[l] < 1):
						rp_val[l] = 1
					#end
				#end
				
				# 
				best_g_old[:] = best_g_val[:]
				
			#end
			
		#end
		
	#end
	
	
	# Print Results
	if (prtoutiter != 0):
		
		# Output to screen
		print("="*80 + "\n")
		print("RANDOM SEED VALUE: %.8f\n" %(rseed))
		print("NUMBER OF ITERATIONS: %d\n" %(k_out))
		print("NUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" %(nfevals))
		print("OBJECTIVE FUNCTION VALUE:")
		print("\tF = %g\n" %(best_f_val))
		if (constraints > 0):
			# Equality Constraints
			print("EQUALITY CONSTRAINTS VALUES:")
			for l in xrange(neqcons):
				print("\tG(%d) = %g" %(l,best_g_val[l]))
			#end
			# Inequality Constraints
			print("\nINEQUALITY CONSTRAINTS VALUES:")
			for l in xrange(neqcons,constraints):
				print("\tH(%d) = %g" %(l,best_g_val[l]))
			#end
		#end
		print("\nLAGRANGIAN MULTIPLIERS VALUES:")
		for l in xrange(constraints):
			print("\tL(%d) = %g" %(l,float(lambda_val[l])))
		#end
		
		print("\nDESIGN VARIABLES VALUES:")
		if (scale == 1):
			x_tmp = (best_x_val[:] * space_halflen) + space_centre
		else:
			x_tmp = best_x_val[:]
		#end
		for m in discrete_i:
			x_tmp[m] = floor(x_tmp[m]+0.5)
		#end
		text = ''
		for j in xrange(dimensions):
			text += ("\tP(%d) = %9.3e\t" %(j,x_tmp[j]))
			if (numpy.mod(j+1,3) == 0):
				text +=("\n")
			#end
		#end
		print text
		print("="*80 + "\n")
	#end
	if (fileout == 1):
		# Output to filename
		ofile.write("\n" + "="*80 + "\n")
		ofile.write("RANDOM SEED VALUE: %.8f\n" %(rseed))
		ofile.write("\nNUMBER OF ITERATIONS: %d\n" %(k_out))
		ofile.write("\nNUMBER OF OBJECTIVE FUNCTION EVALUATIONS: %d\n" %(nfevals))
		ofile.write("\nOBJECTIVE FUNCTION VALUE:\n")
		ofile.write("\tF = %g\n" %(best_f_val))
		if (constraints > 0):
			# Equality Constraints
			ofile.write("\nEQUALITY CONSTRAINTS VALUES:\n")
			for l in xrange(neqcons):
				ofile.write("\tG(%d) = %g\n" %(l,best_g_val[l]))
			#end
			# Inequality Constraints
			ofile.write("\nINEQUALITY CONSTRAINTS VALUES:\n")
			for l in xrange(neqcons,constraints):
				ofile.write("\tH(%d) = %g\n" %(l,best_g_val[l]))
			#end
		#end
		ofile.write("\nLAGRANGIAN MULTIPLIERS VALUES:\n")
		for l in xrange(constraints):
			ofile.write("\tL(%d) = %g\n" %(l,float(lambda_val[l])))
		#end
		
		ofile.write("\nDESIGN VARIABLES VALUES:\n")
		if (scale == 1):
			x_tmp = (best_x_val[:] * space_halflen) + space_centre
		else:
			x_tmp = best_x_val[:]
		#end
		for m in discrete_i:
			x_tmp[m] = floor(x_tmp[m]+0.5)
		#end
		text = ''
		for j in xrange(dimensions):
			text += ("\tP(%d) = %9.3e\t" %(j,x_tmp[j]))
			if (numpy.mod(j+1,3) == 0):
				text +=("\n")
			#end
		#end
		ofile.write(text)
		ofile.write("\n" + "="*80 + "\n")
		
		ofile.close()
	#end
	
	
	# Results
	if (scale == 1):
		opt_x = (best_x_val * space_halflen) + space_centre
	else:
		opt_x = best_x_val
	#end
	for m in discrete_i:
		opt_x[m] = int(floor(opt_x[m] + 0.5))
	#end
	opt_f = best_f_val
	opt_g = best_g_val
	opt_lambda = lambda_val[:]
	
	return opt_x,opt_f,opt_g,opt_lambda,nfevals,'%.8f' %(rseed)
	


# =============================================================================
# chso Function
# =============================================================================
def chso(ND,nc,nec,xtype,x0,lb,ub,bw,HMS,HMCR,PAR,maxIter,printout,rseed,objfunc):
	
	'''
	CHSO function - Python Version of the Constrained Harmony Search Optimizer
	
	Documentation last updated:  October. 16, 2008 - Ruben Perez
	'''
	
	# Set random number seed
	rand = random.Random()
	if rseed == {}:	
		rseed = time.time()
	#end
	
	
	# Initialize
	HM = numpy.zeros((HMS,ND+1), float)
	for i in xrange(HMS):
		for j in xrange(ND):
			HM[i,j] = lb[j] + rand.random()*(ub[j] - lb[j])
		#end
		[f0,gs0] = objfunc(HM[i,:-1])
		HM[i,ND] = f0
	#end
	
	# Print Initial Header
	if (printout == 1):
		#print(' Iteration   Func-count     min f(x)')
		print(' Iteration   min f(x)');
	#end
	
	
	# Iterations Loop
	x = numpy.zeros(ND,float)
	numFunEvals = 0
	k = 0
	status = 0
	while status != 1:
		
		# New Harmony Improvisation
		for j in xrange(ND):
			
			# 
			if (rand.random() >= HMCR):
				
				# Random Searching
				x[j] = lb[j] + rand.random()*(ub[j] - lb[j])
				
			else:
				
				# Harmony Memory Considering
				x[j] = HM[int(HMS*rand.random()),j]
				
				# Pitch Adjusting
				if (rand.random() <= PAR):
					if (rand.random() > 0.5):
						x[j] = x[j] + rand.random()*((ub[j] - lb[j])/bw[j])
					else:
						x[j] = x[j] - rand.random()*((ub[j] - lb[j])/bw[j])
					#end
				#end
				
			#end
		#end
		
		# 
		[fval,gvals] = objfunc(x)
		numFunEvals += 1
		
		# 
		if (sum(gvals) <= 0):
			
			# Harmony Memory Update
			hmax_num = 0
			hmax = HM[0,ND]
			for i in xrange(HMS):
				if (HM[i,ND] > hmax):
					hmax_num = i
					hmax = HM[i,ND]
				#end
			#end
			
			if (fval < hmax):
				for j in xrange(ND):
					HM[hmax_num,j] = x[j]
				#end
				HM[hmax_num,ND] = fval
			#end
			
			hmin_num = 0
			hmin = HM[0,ND]
			for i in xrange(HMS):
				if (HM[i,ND] < hmin):
					hmin_num = i
					hmin = HM[i,ND]
				#end
			#end
			
			# Print
			if (fval == hmin):
				opt_x = x
				opt_f = fval
				opt_g = gvals
				if (printout == 1):
					#print('%f,%f,%f,%f' %(k,x,fval,numpy.var(numpy.corrcoef(HM).T)))
					print('%i,%f' %(k,fval))
				#end
			#end
			
		#end
		
		# Test Convergence
		if k == maxIter-1:
			if (printout == 1):
				print '\nMaximum number of iterations exceeded\n'
				print 'increase OPTIONS.MaxIter\n'
			#end
			status = 1
		else:
			k += 1
		#end
		
	#end
	
	# Print
	if (printout == 1):
		print '\nNumber of function evaluations = %f\n' %(numFunEvals)
	#end
	
	return opt_x,opt_f,opt_g,numFunEvals,'%.8f' %(rseed)
	

#==============================================================================
# Optimizers Test
#==============================================================================
if __name__ == '__main__':
	
	print 'Testing ...'
	
	# Test alpso
	alhso = alhso()
	print alhso
	
