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

    # Operations with scalars

    def __mul__(self, other):
        if isinstance(other, (float, complex, int)):
            return FermiOp(self.orbs, self.pFrom, self.pTo, data = other*self.data)
        else:
            raise TypeError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __rmul__(self, other):
        if isinstance(other, (float, complex, int)):
            return FermiOp(self.orbs, self.pFrom, self.pTo, data = other*self.data)
        else:
            raise TypeError("Argument for multiplication must be numeric. For operator multiplication, use the matmul operator '@' instead.")

    def __truediv__(self, other):
        if isinstance(other, (float, complex, int)):
            return FermiOp(self.orbs,  self.pFrom, self.pTo, data = self.data / other)
        else:
            raise ValueError("Argument for division must be numeric.")

    # Operations with other FermiOp

    def __add__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.pFrom == other.pFrom) & (self.pTo == other.pTo))
            new_op = FermiOp(self.orbs, self.pFrom, self.pTo)
            new_op.data = self.data + other.data
            return new_op
        else:
            raise TypeError("Addition implemented only for objects of same type")

    def __sub__(self, other):
        if type(self) == type(other):
            assert ((self.orbs == other.orbs) & (self.pFrom == other.pFrom) & (self.pTo == other.pTo))
            new_op = FermiOp(self.orbs, self.pFrom, self.pTo)
            new_op.data = self.data - other.data
            return new_op
        else:
            raise TypeError("Addition implemented only for objects of same type")
    
    def __matmul__(self, other):
        if type(self) == type(other):
            assert((self.orbs == other.orbs) and (self.shape[1] == other.shape[0]))
            mul_data = self.data @ other.data
            return FermiOp(self.orbs, other.pFrom, self.pTo, data = mul_data)
        else:
            return NotImplemented 
            #raise TypeError("Operator multiplication only implemented for:\n\t\t * FermiOp @ FermiOp --> FermiOp\n\t\t * FermiState @ FermiOp --> FermiOp\n\t\t * FermiOp @ FermiState --> Complex Number ")
    
    # Hermitian conjugate:
    def dagger(self):
        return FermiOp(self.orbs, self.pTo, self.pFrom, self.data.conj().T)
    
    @property
    def T(self):
        return self.dagger()
