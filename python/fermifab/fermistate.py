import numpy as np
from scipy.special import binom
from fermifab.fermiop import *

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

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        state_info = "Fermi State (orbs == {0}, N == {1})".format(self.orbs,self.N)
        data_info = self.data.__repr__()#.replace("array(","      ").replace(")"," ")
        return state_info + "\n\nVector representation w.r.t. ordered Slater basis:\n\n" + data_info


    def __add__(self, other):
        assert ((self.orbs == other.orbs) & (self.N == other.N))
        new_state = FermiState(self.orbs, self.N)
        new_state.data = self.data + other.data
        return new_state

    def __sub__(self, other):
        assert ((self.orbs == other.orbs) & (self.N == other.N))
        new_state = FermiState(self.orbs, self.N)
        new_state.data = self.data - other.data
        return new_state

    def dagger(self):
        return FermiOp(self.orbs, self.N, 0, self.data.conj().reshape((1, -1)))
    
    @property
    def T(self):
        return self.dagger()

    def __mul__(self, other):
        if type(other) in [float, int, complex]:
            return FermiState(self.orbs, self.N, data = other*self.data)
        elif type(other) == FermiOp:
            assert other.shape[0] == 1
            mul_data = self.data.reshape((-1, 1)) @ other.data
            return FermiOp(self.orbs, other.pFrom, self.N, data = mul_data)
        else: raise

    def __rmul__(self, other):
        if type(other) in [float, int, complex]:
            return FermiState(self.orbs, self.N, data = other*self.data)
            
        # We implement the FermiOp * FermiState operation here
        # to avoid a circular import.
        elif type(other) == FermiOp:
            assert((self.orbs == other.orbs) and (other.shape[1] == len(self)))
            mul_data = other.data @ self.data
            if len(mul_data) == 1:
                return complex(mul_data[0])
            else:
                return FermiState(self.orbs, other.pTo, data = mul_data)
        else: raise