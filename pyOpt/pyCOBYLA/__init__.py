try:
    from .pyCOBYLA import COBYLA
    __all__ = ['COBYLA']
except:
    __all__ = []
