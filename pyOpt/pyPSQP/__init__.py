#!/usr/local/bin/python

try:
    from pyPSQP import PSQP
    __all__ = ['PSQP']
except:
    __all__ = []
#end
