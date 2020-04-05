import numpy as np
from scipy.special import binom
from .fermiop import FermiOp

__all__ = ['FermiState']


class FermiState(object):

    def __init__(self, orbs, N, data=None):
        """
        Construct a Fermi state.

        Args:
            orbs:  number of orbitals
            N:     number of particles
            data:  state vector (optional)
        """
        # TODO:
        # Support for lists of orbs and N
        # if not hasattr(orbs, '__len__'): (orbs,)
        # if not hasattr(N, '__len__'): (N,)
        # self.orbs = np.array(orbs)
        # self.N = np.array(N)
        self.orbs = orbs
        self.N = N
        assert np.all((self.orbs >= self.N) and (self.N >= 0))
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
        if isinstance(other, (float, complex, int)):
            return FermiState(self.orbs, self.N, data=(other*self.data))
        else:
            raise ValueError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __rmul__(self, other):
        if isinstance(other, (float, complex, int)):
            return FermiState(self.orbs, self.N, data=(other*self.data))
        else:
            raise ValueError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __truediv__(self, other):
        if isinstance(other, (float, complex, int)):
            return FermiState(self.orbs, self.N, data=(self.data / other))
        else:
            raise ValueError("Argument for division must be numeric.")

    # Operations with other fermi states
    def __add__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.N == other.N))
            return FermiState(self.orbs, self.N, data=(self.data + other.data))
        else:
            raise TypeError("Addition implemented only for objects of same type")

    def __sub__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.N == other.N))
            return FermiState(self.orbs, self.N, data=(self.data - other.data))
        else:
            raise TypeError("Subtraction implemented only for objects of same type")

    # Operations with FermiOps
    def __matmul__(self, other):
        if type(other) == FermiOp:
            assert other.shape[0] == 1
            # construct outer Kronecker product
            return FermiOp(self.orbs, other.pFrom, self.N, data=(self.data.reshape((-1, 1)) @ other.data))
        else:
            raise TypeError("Operator multiplication only implemented for:\n\t\t * FermiOp @ FermiOp --> FermiOp\n\t\t * FermiState @ FermiOp --> FermiOp\n\t\t * FermiOp @ FermiState --> Complex Number ")

    def __rmatmul__(self, other):
        # We implement the FermiOp * FermiState operation here
        # to avoid a circular import.
        if type(other) == FermiOp:
            assert self.orbs == other.orbs
            return FermiState(self.orbs, other.pTo, data=(other.data @ self.data))
        else:
            raise TypeError("Operator multiplication only implemented for:\n\t\t * FermiOp @ FermiOp --> FermiOp\n\t\t * FermiState @ FermiOp --> FermiOp\n\t\t * FermiOp @ FermiState --> Complex Number ")

    # Hermitian conjugate (construct 'bra' as operator)

    def dagger(self):
        return FermiOp(self.orbs, self.N, 0, self.data.conj().reshape((1, -1)))

    @property
    def H(self):
        return self.dagger()

    def copy(self):
        return FermiState(self.orbs, self.N, self.data.copy())

