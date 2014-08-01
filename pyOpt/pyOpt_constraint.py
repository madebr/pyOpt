#!/usr/bin/env python
'''
pyOpt_constraint

Holds the Python Design Optimization Classes (base and inherited).

Copyright (c) 2008-2014 by pyOpt Developers
All rights reserved.
Revision: 1.1   $Date: 08/05/2008 21:00$


Developers:
-----------
- Dr. Ruben E. Perez (RP)
- Mr. Peter W. Jansen (PJ)

History
-------
    v. 1.0  - Initial Class Creation (RP, 2008)
    v. 1.1  - Pretty Print of Optimization Problems (PJ, 2008)
'''

__version__ = '$Revision: $'

'''
To Do:
    - 
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os, sys
import pdb

# =============================================================================
# External Python modules
# =============================================================================
#import external

# =============================================================================
# Extension modules
# =============================================================================
#import extension

# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity


# =============================================================================
# Constraint Class
# =============================================================================
class Constraint(object):
    
    '''
    Optimization Constraint Class
    '''
    
    def __init__(self, name, type='i', *args, **kwargs):
        
        '''
        Constraint Class Initialization
        
        **Arguments:**
        
        - name -> STR: Variable Name
        
        **Keyword arguments:**
        
        - type -> STR: Variable Type ('i'-inequality, 'e'-equality), *Default* = 'i'
        - lower -> INT: Variable Lower Value
        - upper -> INT: Variable Upper Value
        - choices -> DICT: Variable Choices
        
        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        '''
        
        # 
        self.name = name
        self.type = type[0].lower()
        self.value = 0.0
        if (type[0].lower() == 'i'):
            self.upper = 0.0	#float(inf) 
            self.lower = -float(inf)
            for key in kwargs.keys():
                if (key == 'lower'):
                    self.lower = float(kwargs['lower'])
                #else:
                    #self.lower = -float(inf)
                #end
                if (key == 'upper'):
                    self.upper = float(kwargs['upper'])
                #else:
                    #self.upper = float(inf) 
                #end
            #end
        elif (type[0].lower() == 'e'):
            if (kwargs.has_key('equal')):
                self.equal = float(kwargs['equal'])
            else:
                self.equal = 0.0
            #end
        else:
            raise IOError('Constraint type not understood -- use either i(nequality) or e(quality)')
        #end
        
        #if (kwargs['nvars']):
        #	self.sensitivity = numpy.zeros(kwargs['nvars'],float)
        ##end
        
        
    def ListAttributes(self):
        
        '''
        Print Structured Attributes List
        
        Documentation last updated:  March. 10, 2008 - Ruben E. Perez
        '''
        
        ListAttributes(self)
        
        
    def __str__(self):
        
        '''
        Print Constraint
        
        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        '''
        
        if (self.type == 'e'):
            return ( '	    Name        Type'+' '*25+'Bound\n'+'	 '+str(self.name).center(9) +'    e %23f = %5.2e\n' %(self.value,self.equal))
        if (self.type == 'i'):
            return ( '	    Name        Type'+' '*25+'Bound\n'+'	 '+str(self.name).center(9) +'	  i %15.2e <= %8f <= %8.2e\n' %(self.lower,self.value,self.upper))
    


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
# Constraint Test
#==============================================================================
if __name__ == '__main__':
    
    print 'Testing ...'
    
    # Test Constraint
    con = Constraint('g')
    con.ListAttributes()
    
