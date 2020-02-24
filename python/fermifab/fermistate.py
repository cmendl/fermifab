import numpy as np
from scipy.special import binom
from fermiop import *

__all__ = ['FermiState']


class FermiState(object):

    def __init__(self, orbs, N, data=None):
        self.orbs = orbs
        self.N = N
        if data is not None:
            # TODO: support lists of 'orbs' and 'N'
            assert len(data) == int(binom(orbs, N))
            self.data = np.array(data)
        else:
            self.data = np.zeros(int(binom(orbs, N)), dtype=complex)

    def __add__(self,other):
        assert ((self.orbs == other.orbs) & (self.N == other.N))
        new_state = FermiState(self.orbs,self.N)
        new_state.data = self.data + other.data
        new_state.data /= np.linalg.norm(new_state.data)
        return new_state

    def dagger(self):
        return FermiOp(self.orbs,self.N,0,self.data.conj().reshape((1,-1)))
    
    @property
    def T(self):
        return self.dagger()

    def __mul__(self,other):
        assert self.orbs == other.orbs
        if type(other) == FermiOp:
            mul_data = self.data.reshape((-1,1)) @ other.data
            return FermiOp(self.orbs, other.pFrom, self.N,data = mul_data)
        else: raise