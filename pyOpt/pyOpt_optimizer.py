#!/usr/bin/env python
"""pyOpt_optimizer.

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
            - Bug Fix in setOption (PJ, 2008)
"""

__version__ = '$Revision: $'

'''
To Do:
    - deal with g>=0 & g<=0 depending on optimizer !
    - add timing routine in optimizer class ?? better leave @ specific opt level ?
    - parameters class ?
    - How to handle external gradients calculations ?
    - Implement Restart Capability ?
    - self tests ?
    - add Optimizer Info method
'''

# =============================================================================
# Standard Python modules
# =============================================================================
import os
import sys

# =============================================================================
# Extension modules
# =============================================================================
from .pyOpt_history import History
from .pyOpt_optimization import Optimization

# =============================================================================
# External Python modules
# =============================================================================
#import external


# =============================================================================
# Misc Definitions
# =============================================================================
inf = 10.E+20  # define a value for infinity


# =============================================================================
# Optimizer Class
# =============================================================================
class Optimizer:
    """Abstract Class for Optimizer Object."""

    def __init__(self, name=None, category=None, def_options=None, informs=None, *args, **kwargs):
        """Optimizer Class Initialization.

        **Keyword arguments:**

        - name -> STR: Optimizer name, *Default* = {}
        - category -> STR: Optimizer category, *Default* = {}
        - def_options -> DICT: Deafult options, *Default* = {}
        - informs -> DICT: Calling routine informations texts, *Default* = {}

        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        """

        self.name = name or {}
        self.category = category or {}
        self.options = {}
        self.options['defaults'] = def_options or {}
        self.informs = informs or {}

        # Initialize Options
        def_keys = def_options.keys()
        for key in def_keys:
            self.options[key] = def_options[key]

        koptions = kwargs.pop('options',{})
        kopt_keys = koptions.keys()
        for key in kopt_keys:
            self.setOption(key,koptions[key])

    def __solve__(self, opt_problem, *args, **kwargs):
        """Run Optimizer (Optimizer Specific Routine)

        **Keyword arguments:**

        - opt_problem -> INST: Optimization problem instance, *Default* = {}

        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        """
        raise NotImplementedError()


    def __call__(self, opt_problem, *args, **kwargs):
        """Run Optimizer (Calling Routine)

        **Keyword arguments:**

        - opt_problem -> INST: Optimization problem instance

        Additional arguments and keyword arguments are passed to the objective function call

        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        """

        # Check Optimization Problem
        if not isinstance(opt_problem,Optimization):
            try:
                hasattr(opt_problem,'_constraints')
            except:
                raise ValueError("Input is not a Valid Optimization Problem Instance\n")

        # Check order of Constraints
        last_eq = 0
        first_ieq = -1
        if (len(opt_problem._constraints.keys()) > 0):
            for key in opt_problem._constraints.keys():
                if opt_problem._constraints[key].type == 'e':
                    last_eq = int(key)
                elif opt_problem._constraints[key].type == 'i':
                    if (first_ieq == -1):
                        first_ieq = int(key)

            if (last_eq > first_ieq) and (first_ieq != -1):
                print('WARNING - Equality Constraints should be defined BEFORE Inequality Constraints')

        # Solve Optimization Problem
        return self.__solve__(opt_problem, *args, **kwargs)

    def _on_setOption(self, name, value):
        """Set Optimizer Option Value (Optimizer Specific Routine)

        **Arguments:**

        - name -> STR: Option name
        - value ->   : Option value

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """

        raise NotImplementedError()


    def setOption(self, name, value=None):
        """Set Optimizer Option Value (Calling Routine)

        **Arguments:**

        - name -> STR: Option Name

        **Keyword arguments:**

        - value -> FLOAT/INT/BOOL: Option Value, *Default* = None

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """

        def_options = self.options['defaults']
        if name in def_options:
            if (type(value) == def_options[name][0]):
                self.options[name] = [type(value),value]
            else:
                raise OSError('Incorrect ' + repr(name) + ' value type')
        else:
            raise OSError(repr(name) + ' is not a valid option name')

        self._on_setOption(name, value)

    def _on_getOption(self, name):
        """Get Optimizer Option Value (Optimizer Specific Routine)

        **Arguments:**

        - name -> STR: Option name

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """
        raise NotImplementedError()

    def getOption(self, name):
        """Get Optimizer Option Value (Calling Routine)

        **Arguments:**

        - name -> STR: Option name

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """

        def_options = self.options['defaults']
        if name in def_options:
            return self.options[name][1]
        else:
            raise OSError(repr(name) + ' is not a valid option name')

        self._on_getOption(name)


    def _on_getInform(self, info):
        """Get Optimizer Result Information (Optimizer Specific Routine)

        **Arguments:**

        - info -> STR: Information key

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """
        raise NotImplementedError()

    def getInform(self, infocode=None):
        """Get Optimizer Result Information (Calling Routine)

        **Keyword arguments:**

        - infocode -> INT: information code key

        Documentation last updated:  Feb. 07, 2011 - Peter W. Jansen
        """

        if (infocode == None):
            return self.informs
        else:
            return self._on_getInform(infocode)

    def _on_flushFiles(self):
        """Flush Output Files (Optimizer Specific Routine)

        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """

        raise NotImplementedError()

    def flushFiles(self):
        """Flush Output Files (Calling Routine)

        Documentation last updated:  August. 09, 2009 - Ruben E. Perez
        """

        self._on_flushFiles()

    def _setHistory(self, probname, store_hst, hot_start, def_fname):
        """Setup Optimizer History and/or Hot-start instances.

        **Arguments:**

        - probname  -> STR: Optimization problem name
        - store_hst -> BOOL/STR: Flag/filename to store optimization history
        - hot_start -> BOOL/STR: Flag/filename to read optimization history
        - def_fname -> STR: Default file name

        Documentation last updated:  Oct. 12, 2011 - Peter W. Jansen
        """

        myrank = self.myrank

        hos_file = None
        log_file = None
        tmp_file = False
        if (myrank == 0):
            if isinstance(store_hst,str):
                if isinstance(hot_start,str):
                    if (store_hst == hot_start):
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(store_hst+'_tmp', 'w', self, probname)
                        tmp_file = True
                    else:
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(store_hst, 'w', self, probname)

                    self.sto_hst = True
                    self.hot_start = True
                elif hot_start:
                    hos_file = History(store_hst, 'r', self)
                    log_file = History(store_hst+'_tmp', 'w', self, probname)
                    self.sto_hst = True
                    self.hot_start = True
                    tmp_file = True
                else:
                    log_file = History(store_hst, 'w', self, probname)
                    self.sto_hst = True
                    self.hot_start = False

            elif store_hst:
                if isinstance(hot_start,str):
                    if (hot_start == def_fname):
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(def_fname+'_tmp', 'w', self, probname)
                        tmp_file = True
                    else:
                        hos_file = History(hot_start, 'r', self)
                        log_file = History(def_fname, 'w', self, probname)

                    self.sto_hst = True
                    self.hot_start = True
                elif hot_start:
                    hos_file = History(def_fname, 'r', self)
                    log_file = History(def_fname+'_tmp', 'w', self, probname)
                    self.sto_hst = True
                    self.hot_start = True
                    tmp_file = True
                else:
                    log_file = History(def_fname, 'w', self, probname)
                    self.sto_hst = True
                    self.hot_start = False
            else:
                if isinstance(hot_start,str):
                    hos_file = History(hot_start, 'r', self)
                    self.hot_start = True
                elif hot_start:
                    hos_file = History(def_fname, 'r', self)
                    self.hot_start = True
                else:
                    self.hot_start = False

                self.sto_hst = False
        else:
            self.sto_hst = False
            self.hot_start = False

        return hos_file, log_file, tmp_file

    def ListAttributes(self):
        """Print Structured Attributes List.

        Documentation last updated:  March. 24, 2008 - Ruben E. Perez
        """

        ListAttributes(self)



#==============================================================================
#
#==============================================================================
def ListAttributes(self):
        """Print Structured Attributes List.

        Documentation last updated:  March. 24, 2008 - Ruben E. Perez
        """

        print('\n')
        print('Attributes List of: ' + repr(self.__dict__['name']) + ' - ' + self.__class__.__name__ + ' Instance\n')
        self_keys = self.__dict__.keys()
        for key in sorted(self_keys):
            if key != 'name':
                print(str(key) + ' : ' + repr(self.__dict__[key]))

        print('\n')



#==============================================================================
# Optimizer Test
#==============================================================================
if __name__ == '__main__':

    # Test Optimizer
    print('Testing Optimizer...')
    opt = Optimizer()
    opt.ListAttributes()
