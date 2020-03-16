import numpy as np
from fermifab.fermiop import *
from fermifab.fermistate import *

__all__ = ['crand','norm']

def crand(args):
    return 0.5 - np.random.rand(*args) + 1j*(0.5 - np.random.rand(*args))

def norm(FermiObj):
    if type(FermiObj) in [FermiState,FermiOp]:
        return float(np.linalg.norm(FermiObj.data))
    else:
        return TypeError("Argument must be FermiState or FermiOp")