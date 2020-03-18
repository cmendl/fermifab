import numpy as np
from scipy.special import binom
from fermifab.fermiop import *

__all__ = ['FermiState']

class FermiState(object):
    def __init__(self, orbs, N, data=None):
        """Construct a Fermi state

        Parameters:
        ---------
        orbs:   int
                Number of orbitals.

        N:      int
                Number of particles.

        data:   NumPy array, optional
                State vector.
                
        """
        # Support for lists of orbs and N
        # if not hasattr(orbs, '__len__'): (orbs,)
        # if not hasattr(N, '__len__'): (N,)
        # self.orbs = np.array(orbs)
        # self.N = np.array(N)
        self.orbs = orbs
        self.N = N
        assert np.all((self.orbs >= self.N) & (self.N >= 0))
        if data is not None:
            # TODO: support lists of 'orbs' and 'N'
            assert len(data) == int(np.prod(binom(orbs, N)))
            self.data = np.array(data)
        else:
            self.data = np.zeros(int(np.prod(binom(orbs, N))), dtype=complex)
            self.data[0] = 1.

    def __len__(self):
        return len(self.data)

    def __repr__(self):
        state_info = "Fermi State (orbs == {0}, N == {1})".format(self.orbs,self.N)
        data_info = self.data.__repr__()#.replace("array(","      ").replace(")"," ")
        return state_info + "\n\nVector representation w.r.t. ordered Slater basis:\n\n" + data_info

    # Operations with scalars

    def __mul__(self, other):
        if type(other) in [float, int, complex]:
            return FermiState(self.orbs, self.N, data = other*self.data)
        else:
            raise ValueError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __rmul__(self, other):
        if type(other) in [float, int, complex]:
            return FermiState(self.orbs, self.N, data = other*self.data)
        else:
            raise ValueError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __truediv__(self, other):
        if type(other) in [float, int, complex]:
            return FermiState(self.orbs, self.N, data = self.data / other)
        else:
            raise ValueError("Argument for division must be numeric.")

    # Operations with other Fermistates
    def __add__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.N == other.N))
            new_state = FermiState(self.orbs, self.N)
            new_state.data = self.data + other.data
            return new_state
        else:
            raise TypeError("Addition implemented only for objects of same type")

    def __sub__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.N == other.N))
            new_state = FermiState(self.orbs, self.N)
            new_state.data = self.data - other.data
            return new_state
        else:
            raise TypeError("Addition implemented only for objects of same type")

    # Operations with FermiOps
    def __matmul__(self, other):
        if type(other) == FermiOp:
            assert other.shape[0] == 1
            mul_data = self.data.reshape((-1, 1)) @ other.data
            return FermiOp(self.orbs, other.pFrom, self.N, data = mul_data)
        else: 
            raise TypeError("Operator multiplication only implemented for:\n\t\t * FermiOp @ FermiOp --> FermiOp\n\t\t * FermiState @ FermiOp --> FermiOp\n\t\t * FermiOp @ FermiState --> Complex Number ")
    
    def __rmatmul__(self, other):
        # We implement the FermiOp * FermiState operation here
        # to avoid a circular import.
        if type(other) == FermiOp:
            assert((self.orbs == other.orbs) and (other.shape[1] == len(self)))
            mul_data = other.data @ self.data
            if len(mul_data) == 1:
                return complex(mul_data[0])
            else:
                return FermiState(self.orbs, other.pTo, data = mul_data)
        else: 
            raise TypeError("Operator multiplication only implemented for:\n\t\t * FermiOp @ FermiOp --> FermiOp\n\t\t * FermiState @ FermiOp --> FermiOp\n\t\t * FermiOp @ FermiState --> Complex Number ")
    
    # Hermitian conjugate

    def dagger(self):
        return FermiOp(self.orbs, self.N, 0, self.data.conj().reshape((1, -1)))
    
    @property
    def T(self):
        return self.dagger()

