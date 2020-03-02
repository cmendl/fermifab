import numpy as np
from scipy.special import binom


__all__ = ['FermiOp']

class FermiOp(object):

    def __init__(self, orbs, pFrom, pTo, data=None):
        self.orbs = orbs
        self.pFrom = pFrom
        self.pTo = pTo
        self.shape = (int(binom(orbs, pTo)), int(binom(orbs, pFrom)))
        if data is not None:
            assert data.shape == self.shape
            self.data = np.array(data)
        else:
            self.data = np.zeros(self.shape, dtype=complex)

    def __repr__(self):
        state_info = "Fermi Operator wedge^{0} H -> wedge^{1} H (orbs == {2})".format(self.pFrom, self.pTo, self.orbs)
        data_info = self.data.__repr__()
        return state_info + "\n\nMatrix representation w.r.t. ordered Slater basis:\n\n" + data_info

    def __mul__(self, other):
        if type(other) in [float, int, complex]:
            return FermiOp(self.orbs, self.pFrom, self.pTo, data = other*self.data)
        elif type(other) == type(self):
            assert((self.orbs == other.orbs) and (self.shape[1] == other.shape[0]))
            mul_data = self.data @ other.data
            return FermiOp(self.orbs, other.pFrom, self.pTo, data = mul_data)
        else: 
            return NotImplemented 
    
    def __rmul__(self,other):
        if type(other) in [float, int, complex]:
            return FermiOp(self.orbs, self.pFrom, self.pTo, data = other*self.data)
        else:
            return NotImplemented