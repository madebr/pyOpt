#!/usr/bin/env python
"""pyOpt_constraint.

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
"""

__version__ = '$Revision: $'

inf = 10.E+20  # define a value for infinity


class Constraint:
    """Optimization Constraint Class."""

    def __init__(self, name, type='i', *args, **kwargs):

        """Constraint Class Initialization.

        **Arguments:**

        - name -> STR: Variable Name

        **Keyword arguments:**

        - type -> STR: Variable Type ('i'-inequality, 'e'-equality), *Default* = 'i'
        - lower -> INT: Variable Lower Value
        - upper -> INT: Variable Upper Value
        - choices -> DICT: Variable Choices

        Documentation last updated:  Feb. 03, 2011 - Peter W. Jansen
        """

        #
        self.name = name
        self.type = type[0].lower()
        self.value = 0.0
        if type[0].lower() == 'i':
            self.upper = 0.0  # float(inf)
            self.lower = -float(inf)
            for key in kwargs.keys():
                if key == 'lower':
                    self.lower = float(kwargs['lower'])
                # else:
                # self.lower = -float(inf)
                if key == 'upper':
                    self.upper = float(kwargs['upper'])
                # else:
                # self.upper = float(inf)
        elif type[0].lower() == 'e':
            if 'equal' in kwargs:
                self.equal = float(kwargs['equal'])
            else:
                self.equal = 0.0
        else:
            raise OSError('Constraint type not understood -- use either i(nequality) or e(quality)')

        # if (kwargs['nvars']):
        #	self.sensitivity = numpy.zeros(kwargs['nvars'],float)

    def ListAttributes(self):

        """Print Structured Attributes List.

        Documentation last updated:  March. 10, 2008 - Ruben E. Perez
        """

        ListAttributes(self)

    def __str__(self):

        """Print Constraint.

        Documentation last updated:  April. 30, 2008 - Peter W. Jansen
        """

        if self.type == 'e':
            return ('	    Name        Type' + ' ' * 25 + 'Bound\n' + '	 ' + str(self.name).center(
                9) + f'    e {self.value:23f} = {self.equal:5.2e}\n')
        if self.type == 'i':
            return ('	    Name        Type' + ' ' * 25 + 'Bound\n' + '	 ' + str(self.name).center(
                9) + f'	  i {self.lower:15.2e} <= {self.value:8f} <= {self.upper:8.2e}\n')


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


# ==============================================================================
# Constraint Test
# ==============================================================================
if __name__ == '__main__':
    print('Testing ...')

    # Test Constraint
    con = Constraint('g')
    con.ListAttributes()
