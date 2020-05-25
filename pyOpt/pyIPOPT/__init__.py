try:
    from .pyIPOPT import IPOPT
    __all__ = ['IPOPT']
except:
    __all__ = []
