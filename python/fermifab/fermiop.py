import numpy as np
from scipy.special import binom


__all__ = ['FermiOp']

class FermiOp(object):
    

    def __init__(self, orbs, pFrom, pTo, data=None):
        """Construct a Fermi operator between two Hilbert spaces.

        Parameters:
        ---------
        orbs:   int
                Number of orbitals.

        pFrom:  int
                Number of particles in the input Hilbert space.

        pTo:    int
                Number of particles in the output Hilbert space.

        data:   NumPy array, optional
                Operator matrix.

        """
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

    def __add__(self, other):
        assert ((self.orbs == other.orbs) & (self.pFrom == other.pFrom) & (self.pTo == other.pTo))
        new_state = FermiOp(self.orbs, self.pFrom, self.pTo)
        new_state.data = self.data + other.data
        return new_state

    def __sub__(self, other):
        assert ((self.orbs == other.orbs) & (self.pFrom == other.pFrom) & (self.pTo == other.pTo))
        new_state = FermiOp(self.orbs, self.pFrom, self.pTo)
        new_state.data = self.data - other.data
        return new_state

    def dagger(self):
        return FermiOp(self.orbs, self.pTo, self.pFrom, self.data.conj().T)
    
    @property
    def T(self):
        return self.dagger()

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