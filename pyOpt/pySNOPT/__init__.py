try:
    from .pySNOPT import SNOPT
    __all__ = ['SNOPT']
except:
    __all__ = []
