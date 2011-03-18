#!/usr/local/bin/python

try:
    from pyMIDACO import MIDACO
    __all__ = ['MIDACO']
except:
    __all__ = []
#end
