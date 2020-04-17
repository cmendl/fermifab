from scipy.sparse import csr_matrix
from .fermiop import FermiOp
import fermifab.kernel

__all__ = ['tensor_op']


def tensor_op(op, N):
    """
    Calculate the matrix representation of the N-fold tensor product of an operator.

    Args:
        op: quantum operator of type `FermiOp`, with `pFrom` and `pTo` equal to 1
        N:  number of tensor factors

    Returns:
        FermiOp: N-fold tensor product
    """
    if op.pFrom != 1 or op.pTo != 1:
        raise ValueError('operator particle numbers must be equal to 1')
    dims, val, ind = fermifab.kernel.tensor_op(op.data, N)
    # finally convert to dense matrix, for simplicity
    AN = csr_matrix((val, (ind[:, 0], ind[:, 1])), shape=dims).todense()
    return FermiOp(op.orbs, N, N, data=AN)
