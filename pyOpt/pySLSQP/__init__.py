try:
    from .pySLSQP import SLSQP
    __all__ = ['SLSQP']
except:
    __all__ = []
