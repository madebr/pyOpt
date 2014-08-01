#!/usr/bin/env python
'''
pyOpt_optimization

Holds the Python Design Optimization Classes (base and inherited).

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.4   $Date: 22/06/2009 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
    v. 1.0  - Initial Class Creation (RP, 2008)
            - Added Set Variable and Constraint Groups (PJ, 2008)
    v. 1.1  - Pretty Print of Optimization Problems (PJ, 2008)
    v. 1.2  - Added Solution Class (PJ, 2008)
            - Added File Writing Support (PJ, 2008)
    v. 1.3  - Minor Fixes and Functionality Updates (RP, 2008)
    v. 1.4  - Added Variables Groups Handling (PJ,RP 2009)
'''

__version__ = '$Revision: $'

'''
To Do:
    - add variable group error when groups have the same name
    - add method for addVar2Group
    - pickle wrapping ?!
    - save class __str__ info to file (text/TeX) ?
    - warm start from other opts?
    - class for core sensitivity?
    - class for history?
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb

# =============================================================================
# External Python modules
# =============================================================================
import numpy

# =============================================================================
# Extension modules
# =============================================================================
from pyOpt import Variable
from pyOpt import Objective
from pyOpt import Constraint
from pyOpt import Parameter

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity


# =============================================================================
# Optimization Class
# =============================================================================
class Optimization(object):
    
    '''
    Optimization Problem Class
    '''
    
    def __init__(self, name, obj_fun, var_set=None, obj_set=None, con_set=None, use_groups=False, *args, **kwargs):
        
        '''
        Optimization Problem Class Initialization
        
        **Arguments:**
        
        - name -> STR: Solution name
        - opt_func -> FUNC: Objective function
        
        **Keyword arguments:**
        
        - var_set -> INST: Variable set, *Default* = None
        - obj_set -> INST: Objective set, *Default* = None
        - con_set -> INST: Constraints set, *Default* = None
        - use_groups -> BOOL: Use of group identifiers flag, *Default* = False
        
        Documentation last updated:  May. 23, 2011 - Ruben E. Perez
        '''
        
        # 
        self.name = name
        self.obj_fun = obj_fun
        self.use_groups = use_groups
        
        # Initialize Variable Set
        if var_set is None:
            self._variables = {}
        else:
            self._variables = var_set
        #end
        self._vargroups = {}
        
        # Initialize Objective Set
        if obj_set is None:
            self._objectives = {}
        else:
            self._objectives = obj_set
        #end
        
        # Initialize Constraint Set
        if con_set is None:
            self._constraints = {}
        else:
            self._constraints = con_set
        #end
        
        ## Initialize Parameter Set
        #if par_set is None:
        #    self._parameters = {}
        #else:
        #    self._parameters = par_set
        ##end
        #self._pargroups = {}
        
        # Initialize Solution Set
        self._solutions = {}
        
        
    def getVar(self, i):
        
        '''
        Get Variable *i* from Variables Set
        
        **Arguments:**
        
        - i -> INT: Variable index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Variable index must be an integer >= 0.")
        #end
        
        # 
        return self._variables[i]
        
        
    def addVar(self, *args, **kwargs):
        
        '''
        Add Variable into Variables Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        # 
        id = self.firstavailableindex(self._variables)
        self.setVar(id,*args,**kwargs)
        
        # 
        tmp_group = {}
        tmp_group[self._variables[id].name] = id
        self._vargroups[self.firstavailableindex(self._vargroups)] = {'name':self._variables[id].name,'ids':tmp_group}
        
        
    def addVarGroup(self, name, nvars, type='c', value=0.0, **kwargs):
        
        '''
        Add a Group of Variables into Variables Set
        
        **Arguments:**
        
        - name -> STR: Variable Group Name
        - nvars -> INT: Number of variables in group
        
        **Keyword arguments:**
        
        - type -> STR: Variable type ('c'-continuous, 'i'-integer, 'd'-discrete), *Default* = 'c'
        - value ->INT/FLOAT: Variable starting value, *Default* = 0.0
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        #ngroups = len(self._vargroups)
        #for j in xrange(ngroups):
        #    if (self._vargroups[j]['name'] == name):
        #        raise IOError('Variables group names should be distinct\n')
        #    #end
        ##end
        
        #
        type = [type]*nvars
        
        if isinstance(value,list) or isinstance(value,numpy.ndarray):
            value = value
        elif isinstance(value,int):
            value = [value]*nvars
        elif isinstance(value,float):
            value = [value]*nvars
        else:
            raise IOError('Variable type for value not understood - use float, int or list\n')
        #end
        
        lower = [-inf]*nvars
        upper = [inf]*nvars
        choices = ['']*nvars
        
        for key in kwargs.keys():
            if (key == 'lower'):
                if isinstance(kwargs['lower'],float):
                    lower = [kwargs['lower']]*nvars
                elif isinstance(kwargs['lower'],int):
                    lower = [kwargs['lower']]*nvars
                elif isinstance(kwargs['lower'],(list,numpy.ndarray)):
                    if len(kwargs['lower']) != nvars:
                        for i in xrange(len(kwargs['lower'])):
                            lower[i] = kwargs['lower'][i]
                        #end
                    else:
                        lower = kwargs['lower']
                    #end
                else:
                    raise IOError('Variable type for lower bound not understood - use float, int or list\n')
                #end
            elif (key == 'upper'):
                if isinstance(kwargs['upper'],float):
                    upper = [kwargs['upper']]*nvars
                elif isinstance(kwargs['upper'],int):
                    upper = [kwargs['upper']]*nvars
                elif isinstance(kwargs['upper'],(list,numpy.ndarray)):
                    if len(kwargs['upper']) != nvars:
                        for i in xrange(len(kwargs['upper'])):
                            upper[i] = kwargs['upper'][i]
                        #end
                    else:
                        upper = kwargs['upper']
                    #end
                else:
                    raise IOError('Variable type for upper bound not understood - use float, int or list\n')
                #end
            #end
            if  (key == 'choices'):
                choices = [kwargs['choices']]*nvars
            #end
        #end
        
        tmp_group = {}
        for var in xrange(nvars):
            tmp_name = name +'_%s' %(var)
            id = self.firstavailableindex(self._variables)
            self.setVar(id, tmp_name, type[var], value[var], lower=lower[var], upper=upper[var], choices=choices[var])
            tmp_group[tmp_name] = id
        #end
        self._vargroups[self.firstavailableindex(self._vargroups)] = {'name':name,'ids':tmp_group}
        
        
    def setVar(self, i, *args, **kwargs):
        
        '''
        Set Variable *i* into Variables Set
        
        **Arguments:**
        
        - i -> INT: Variable index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        if (len(args) > 0) and isinstance(args[0], Variable):
            self._variables[i] = args[0]
        else:
            try:
                self._variables[i] = Variable(*args,**kwargs)
            except IOError, (error):
                raise IOError("%s" %(error))
            except:
                raise ValueError("Input is not a Valid for a Variable Object instance\n")
            #end
        #end
        
        
    def delVar(self, i):
        
        '''
        Delete Variable *i* from Variables Set
        
        **Arguments:**
        
        - i -> INT: Variable index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Variable index must be an integer >= 0.")
        #end
        
        # 
        del self._variables[i]
        
        # 
        #ngroups = len(self._vargroups)
        for j in self._vargroups.keys():
            keys = self._vargroups[j]['ids']
            nkeys = len(keys)
            for key in keys:
                if (self._vargroups[j]['ids'][key] == i):
                    del self._vargroups[j]['ids'][key]
                    if (nkeys == 1):
                        del self._vargroups[j]
                    #end
                    return
                #end
            #end
        #end
        
        
    def delVarGroup(self, name):
        
        '''
        Delete Variable Group *name* from Variables Set
        
        **Arguments:**
        
        - name -> STR: Variable group name
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        ngroups = len(self._vargroups)
        for j in xrange(ngroups):
            if (self._vargroups[j]['name'] == name):
                keys = self._vargroups[j]['ids']
                for key in keys:
                    id = self._vargroups[j]['ids'][key]
                    del self._variables[id]
                #end
                del self._vargroups[j]
            #end
        #end
        
        
    def getVarSet(self):
        
        '''
        Get Variables Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        return self._variables
        
        
    def getVarGroups(self):
        
        '''
        Get Variables Groups Set
        
        Documentation last updated:  June. 25, 2009 - Ruben E. Perez
        '''
        
        return self._vargroups
        
        
    def getObj(self, i):
        
        '''
        Get Objective *i* from Objectives Set
        
        **Arguments:**
        
        - i -> INT: Objective index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Objective index must be an integer >= 0.")
        #end
        
        # 
        return self._objectives[i]
        
        
    def addObj(self, *args, **kwargs):
        
        '''
        Add Objective into Objectives Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        # 
        self.setObj(self.firstavailableindex(self._objectives),*args,**kwargs)
        
        
    def setObj(self, i, *args, **kwargs):
        
        '''
        Set Objective *i* into Objectives Set
        
        **Arguments:**
        
        - i -> INT: Objective index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        if (len(args) > 0) and isinstance(args[0], Objective):
            self._objectives[i] = args[0]
        else:
            try:
                self._objectives[i] = Objective(*args,**kwargs)
            except:
                raise ValueError("Input is not a Valid for a Objective Object instance\n")
            #end
        #end
        
        
    def delObj(self, i):
        
        '''
        Delete Objective *i* from Objectives Set
        
        **Arguments:**
        
        - i -> INT: Objective index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Objective index must be an integer >= 0.")
        #end
        
        # 
        del self._objectives[i]
        
        
    def getObjSet(self):
        
        '''
        Get Objectives Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        return self._objectives
        
        
    def getCon(self, i):
        
        '''
        Get Constraint *i* from Constraint Set
        
        **Arguments:**
        
        - i -> INT: Constraint index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Constraint index must be an integer >= 0.")
        #end
        
        # 
        return self._constraints[i]
        
        
    def addCon(self, *args, **kwargs):
        
        '''
        Add Constraint into Constraints Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        # 
        self.setCon(self.firstavailableindex(self._constraints),*args,**kwargs)
        
        
    def addConGroup(self, name, ncons, type='i', **kwargs):
        
        '''
        Add a Group of Constraints into Constraints Set
        
        **Arguments:**
        
        - name -> STR: Constraint group name
        - ncons -> INT: Number of constraints in group
        
        **Keyword arguments:**
        
        - type -> STR: Constraint type ('i'-inequality, 'e'-equality), *Default* = 'i'
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        #
        type_list = [type[0].lower()]*ncons
        
        if (type[0].lower() == 'i'):
            lower = [-inf]*ncons
            upper = [0.0]*ncons			
            for key in kwargs.keys():
                if (key == 'lower'):
                    if isinstance(kwargs['lower'],float):
                        lower = [kwargs['lower']]*ncons
                    elif isinstance(kwargs['lower'],int):
                        lower = [kwargs['lower']]*ncons
                    elif isinstance(kwargs['lower'],(list,numpy.ndarray)):
                        if len(kwargs['lower']) != ncons:
                            for i in xrange(len(kwargs['lower'])):
                                lower[i] = kwargs['lower'][i]
                            #end
                        else:
                            lower = kwargs['lower']
                        #end
                    else:
                        raise IOError('Variable type for lower bound not understood - use float, int or list\n')
                    #end
                elif (key == 'upper'):
                    if isinstance(kwargs['upper'],float):
                        upper = [kwargs['upper']]*ncons
                    elif isinstance(kwargs['upper'],int):
                        upper = [kwargs['upper']]*ncons
                    elif isinstance(kwargs['upper'],(list,numpy.ndarray)):
                        if len(kwargs['upper']) != ncons:
                            for i in xrange(len(kwargs['upper'])):
                                upper[i] = kwargs['upper'][i]
                            #end
                        else:
                            upper = kwargs['upper']
                        #end
                    else:
                        raise IOError('Variable type for upper bound not understood - use float, int or list\n')
                    #end
                #end
            #end
            for con in xrange(ncons):
                tmp_name = name +'_%s' %(con)
                self.setCon(self.firstavailableindex(self._constraints),tmp_name, type_list[con], lower=lower[con], upper=upper[con])
            #end
        elif (type[0].lower() == 'e'):
            equal = [0.0]*ncons
            for key in kwargs.keys():
                if (key == 'equal'):
                    if isinstance(kwargs['equal'],float):
                        equal = [kwargs['equal']]*ncons
                    elif isinstance(kwargs['equal'],int):
                        equal = [kwargs['equal']]*ncons
                    elif isinstance(kwargs['equal'],(list,numpy.ndarray)):
                        if len(kwargs['equal']) != ncons:
                            for i in xrange(len(kwargs['equal'])):
                                lower[i] = kwargs['equal'][i]
                            #end
                        else:
                            equal = kwargs['equal']
                        #end
                    else:
                        raise IOError('Variable type for lower bound not understood - use float, int or list\n')
                    #end
                #end
            #end
            for con in xrange(ncons):
                tmp_name = name +'_%s' %(con)
                self.setCon(self.firstavailableindex(self._constraints),tmp_name, type_list[con], equal=equal[con])
            #end
        #end
        
        
    def setCon(self, i, *args, **kwargs):
        
        '''
        Set Constraint *i* into Constraints Set
        
        **Arguments:**
        
        - i -> INT: Constraint index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        if (len(args) > 0) and isinstance(args[0], Constraint):
            self._constraints[i] = args[0]
        else:
            try:
                self._constraints[i] = Constraint(*args,**kwargs)
            except IOError, (error):
                raise IOError("%s" %(error))
            except:
                raise ValueError("Input is not a Valid for a Constraint Object instance\n")
            #end
        #end
        
        
    def delCon(self, i):
        
        '''
        Delete Constraint *i* from Constraints Set
        
        **Arguments:**
        
        - i -> INT: Constraint index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Constraint index must be an integer >= 0.")
        #end
        
        # 
        del self._constraints[i]
        
        
    def getConSet(self):
        
        '''
        Get Constraints Set
        
        Documentation last updated:  March. 27, 2008 - Ruben E. Perez
        '''
        
        return self._constraints
        
        
    def getSol(self, i):
        
        '''
        Get Solution *i* from Solution Set
        
        **Arguments:**
        
        - i -> INT: Solution index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Solution index must be an integer >= 0.")
        #end
        
        # 
        return self._solutions[i]
        
        
    def addSol(self, *args, **kwargs):
        
        '''
        Add Solution into Solution Set
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        # 
        self.setSol(self.firstavailableindex(self._solutions),*args,**kwargs)
        
        
    def setSol(self,i, *args, **kwargs):
        
        '''
        Set Solution *i* into Solution Set
        
        **Arguments:**
        
        - i -> INT: Solution index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        if (len(args) > 0) and isinstance(args[0], Solution):
            self._solutions[i] = args[0]
        else:
            #try:
            self._solutions[i] = Solution(*args,**kwargs)
            #except:
            #	print args
            #	print kwargs
            #	raise ValueError("Input is not a Valid for a Solution Object instance\n")
            #end
        #end
        
        
    def delSol(self, i):
        
        '''
        Delete *i* Solution from Solutions Set
        
        **Arguments:**
        
        - i -> INT: Solution index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Solution index must be an integer >= 0.")
        #end
        
        # 
        del self._solutions[i]
        
        
    def getSolSet(self):
        
        '''
        Get Solutions Set
        
        Documentation last updated:  May. 07, 2008 - Ruben E. Perez
        '''
        
        return self._solutions
        
        
#    def getPar(self, i):
#        
#        '''
#        Get Parameter *i* from Parameters Set
#        
#        **Arguments:**
#        
#        - i -> INT: Solution index
#        
#        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
#        '''
#        
#        # Check Index
#        if not (isinstance(i,int) and i >= 0):
#            raise ValueError("Parameter index must be an integer >= 0.")
#        #end
#        
#        # 
#        return self._parameters[i]
#        
#        
#    def addPar(self, *args, **kwargs):
#        
#        '''
#        Add Parameter into Parameters Set
#        
#        Documentation last updated:  May. 23, 2008 - Ruben E. Perez
#        '''
#        
#        # 
#        id = self.firstavailableindex(self._parameters)
#        self.setPar(id,*args,**kwargs)
#        
#        # 
#        tmp_group = {}
#        tmp_group[self._parameters[id].name] = id
#        self._pargroups[self.firstavailableindex(self._pargroups)] = {'name':self._parameters[id].name,'ids':tmp_group}
#        
#        
#    def setPar(self, i, *args, **kwargs):
#        
#        '''
#        Set Parameter *i* into Parameters Set
#        
#        **Arguments:**
#        
#        - i -> INT: Parameter index
#        
#        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
#        '''
#        
#        # 
#        if (len(args) > 0) and isinstance(args[0], Parameter):
#            self._parameters[i] = args[0]
#        else:
#            try:
#                self._parameters[i] = Parameter(*args,**kwargs)
#            except IOError, (error):
#                raise IOError("%s" %(error))
#            except:
#                raise ValueError("Input is not a Valid for a Parameter Object instance\n")
#            #end
#        #end
#        
#        
#    def delPar(self, i):
#        
#        '''
#        Delete Parameter *i* from Parameters Set
#        
#        **Arguments:**
#        
#        - i -> INT: Parameter index
#        
#        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
#        '''
#        
#        # Check Index
#        if not (isinstance(i,int) and i >= 0):
#            raise ValueError("Parameter index must be an integer >= 0.")
#        #end
#        
#        # 
#        del self._parameters[i]
#        
#        # 
#        for j in self._pargroups.keys():
#            keys = self._pargroups[j]['ids']
#            nkeys = len(keys)
#            for key in keys:
#                if (self._pargroups[j]['ids'][key] == i):
#                    del self._pargroups[j]['ids'][key]
#                    if (nkeys == 1):
#                        del self._pargroups[j]
#                    #end
#                    return
#                #end
#            #end
#        #end
#        
#        
#    def delParGroup(self, name):
#        
#        '''
#        Delete Parameter Group *name* from Parameters Set
#        
#        **Arguments:**
#        
#        - name -> STR: Parameter group name
#        
#        Documentation last updated:  Sep. 07, 2013 - Ruben E. Perez
#        '''
#        
#        # 
#        ngroups = len(self._pargroups)
#        for j in xrange(ngroups):
#            if (self._pargroups[j]['name'] == name):
#                keys = self._pargroups[j]['ids']
#                for key in keys:
#                    id = self._pargroups[j]['ids'][key]
#                    del self._parameters[id]
#                #end
#                del self._pargroups[j]
#            #end
#        #end
#        
#        
#    def getParSet(self):
#        
#        '''
#        Get Parameter Set
#        
#        Documentation last updated:  May. 23, 2008 - Ruben E. Perez
#        '''
#        
#        return self._parameters
#        
#        
#    def getParGroups(self):
#        
#        '''
#        Get Parameters Groups Set
#        
#        Documentation last updated:  June. 25, 2009 - Ruben E. Perez
#        '''
#        
#        return self._pargroups
        
        
    def firstavailableindex(self, set):
        
        '''
        List First Unused Index from Variable Objects List
        
        **Arguments:**
        
        - set -> LIST: Set to find frist available index of
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # 
        i = 0
        while i in set: 
            i += 1
        #end
        
        return i
        
        
    def ListAttributes(self):
        
        '''
        Print Structured Attributes List
        
        Documentation last updated:  March. 24, 2008 - Ruben E. Perez
        '''
        
        ListAttributes(self)
        
        
    def __str__(self):
        
        '''
        Print Structured Optimization Problem
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        text = '''\nOptimization Problem -- %s\n%s\n
        Objective Function: %s\n\n    Objectives:
        Name        Value        Optimum\n''' %(self.name,'='*80,self.obj_fun.__name__)
        for obj in self._objectives.keys():
            lines = str(self._objectives[obj]).split('\n')
            text += lines[1] + '\n'
        #end
        text += '''\n	Variables (c - continuous, i - integer, d - discrete):
        Name    Type       Value       Lower Bound  Upper Bound\n'''
        for var in self._variables.keys():
            lines = str(self._variables[var]).split('\n')
            text+= lines[1] + '\n'
        #end
        if len(self._constraints.keys()) > 0:
            text += '''\n	Constraints (i - inequality, e - equality):
        Name    Type                    Bounds\n'''
            for con in self._constraints.keys():
                lines = str(self._constraints[con]).split('\n')
                text+= lines[1] + '\n'
            #end
        #end
        
        return (text)
        
        
    def write2file(self, outfile='', disp_sols=False, **kwargs):
        
        '''
        Write Structured Optimization Problem to file
        
        **Keyword arguments:**
        
        - outfile   ->  STR/INST: File name or file instance, *Default* = ''
        - disp_sols ->  BOOL: Display solutions flag, *Default* = False.
        - solutions ->  LIST: List of solution indexes.
        
        Documentation last updated:  May. 9, 2008 - Peter W. Jansen
        '''
        
        # 
        if isinstance(outfile,str):
            if (outfile == ''):
                findir = os.listdir(os.curdir)
                tmpname = self.name.lower()
                tmpname = tmpname.split(' ')
                tmpname = tmpname[0]
                i = 0
                while (tmpname+'.txt') in findir:
                    tmpname = tmpname.rstrip('_%d' %(i-1))
                    tmpname = tmpname + '_' +str(i)
                    i += 1
                #end
                tmpname += '.txt'
                outfile = open(tmpname,'w')
            else:
                outfile = open(outfile,'w')
        elif (not isinstance(outfile,str)) and (not isinstance(outfile,file)):
            raise IOError(repr(outfile) + 'is not a file or filename')
        #end
        ftext = self.__str__()
        outfile.write(ftext)
        if disp_sols or kwargs.has_key('solutions'):
            if kwargs.has_key('solutions'):
                sol_indices = kwargs['solutions']
            else:
                sol_indices = self._solutions.keys()
            #end
            for key in sol_indices:
                soltext = '\n' + self._solutions[key].__str__() 
                outfile.write(soltext)
            #end
        #end
        print 'Data written to file ', outfile.name
        outfile.close()
        
        
    def solution(self, i):
        
        '''
        Get Solution from Solution Set
        
        **Arguments:**
        
        - i -> INT: Solution index
        
        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        '''
        
        # Check Index
        if not (isinstance(i,int) and i >= 0):
            raise ValueError("Solution index must be an integer >= 0.")
        #end
        
        # 
        return self._solutions[i]
    


# =============================================================================
# Solution Class
# =============================================================================
class Solution(Optimization):
    
    '''
    Optimization Solution Class
    '''
    
    def __init__(self, optimizer, name, obj_fun, opt_time, opt_evals, opt_inform, var_set=None, obj_set=None, con_set=None, options_set=None, myrank=0,*args, **kwargs):
        
        '''
        Solution Class Initialization
        
        **Arguments:**
        
        - optimizer -> STR: Optimizer name
        - name -> STR: Optimization problem name
        - opt_time -> FLOAT: Solution total time
        - opt_evals -> INT: Number of function evaluations
        
        **Keyword arguments:**
        
        - var_set -> INST: Variable set, *Default* = {}
        - obj_set -> INST: Objective set, *Default* = {}
        - con_set -> INST: Constraints set, *Default* = {}
        - options_set -> Options used for solution, *Default* = {}
        - myrank -> INT: Process identification for MPI evaluations, *Default* = 0
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # 
        Optimization.__init__(self, name, obj_fun, var_set, obj_set, con_set, *args, **kwargs)
        self.optimizer = optimizer
        self.opt_time = opt_time
        self.opt_evals = opt_evals
        self.opt_inform = opt_inform
        self.options_set = options_set
        self.myrank = myrank
        
        if kwargs.has_key('display_opts'):
            self.display_opt = kwargs['display_opts']
            del kwargs['display_opts']
        else:
            self.display_opt = False
        #end
        self.parameters = kwargs
        
        
    def __str__(self):
        
        '''
        Print Structured Solution
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        text0 = Optimization.__str__(self)
        text1 = ''
        lines = text0.split('\n')
        lines[1] = lines[1][len('Optimization Problem -- '):]
        for i in xrange(5):
            text1 += lines[i] + '\n'
        #end
        if self.display_opt:
            text1 += '\n	Options:\n '
            opt_keys = self.options_set.keys()
            opt_keys.sort()
            for key in opt_keys:
                ns = 25-len(key)
                text1 += '		'+ key +':' + str(self.options_set[key][1]).rjust(ns,'.') + '\n'
            #end
        #end
        text1 += '\n    Solution: \n'
        text1 += ('-'*80) + '\n'
        text1 += '    Total Time: %25.4f\n' %(self.opt_time)
        text1 += '    Total Function Evaluations: %9.0i\n' %(self.opt_evals)
        for key in self.parameters.keys():
            if (isinstance(self.parameters[key],(dict,list,tuple))) and (len(self.parameters[key]) == 0):
                continue
            elif (isinstance(self.parameters[key],numpy.ndarray)) and (0 in (self.parameters[key]).shape):
                continue
            else:
                text1 += '    '+ key +': ' + str(self.parameters[key]).rjust(9) + '\n'
            #end
        #end
        for i in xrange(5,len(lines)):
            text1 += lines[i] + '\n'
        #end
        text1 += ('-'*80) + '\n'
        
        if (self.myrank == 0):
            return text1
        else:
            return ''
        #end
        
        
    def write2file(self, outfile):
        
        '''
        Write Structured Solution to file
        
        **Arguments:**
        
        - outfile -> STR: Output file name
        
        Documentation last updated:  May. 9, 2008 - Peter W. Jansen
        '''
        
        Optimization.write2file(self,outfile,False)
    


#==============================================================================
# 
#==============================================================================
def ListAttributes(self):
    
    '''
    Print Structured Attributes List
    
    Documentation last updated:  March. 24, 2008 - Ruben E. Perez
    '''
    
    print '\n'
    print 'Attributes List of: ' + repr(self.__dict__['name']) + ' - ' + self.__class__.__name__ + ' Instance\n'
    self_keys = self.__dict__.keys()
    self_keys.sort()
    for key in self_keys:
        if key != 'name':
            print str(key) + ' : ' + repr(self.__dict__[key])
        #end
    #end
    print '\n'
    


#==============================================================================
# Optimization Test
#==============================================================================
if __name__ == '__main__':
    
    print 'Testing Optimization...'
    optprob = Optimization('Optimization Problem',{})
    optprob.ListAttributes()
    
