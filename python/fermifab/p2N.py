import numpy as np
from .fermiop import FermiOp
from .rdm import construct_rdm_kernel

__all__ = ['p2N']


def p2N(h, N):
    """
    Calculate N-body from p-body operator:
    H = sum_ij h_ij a^dagger_i a_j

    Args:
        h: p-body operator
        N: target particle number

    Returns:
        FermiOp: N-body operator generated from h
    """
    assert type(h) == FermiOp
    assert not (np.array(h.pFrom) - np.array(h.pTo)).any()

    # kernel matrices for current configuration
    K = construct_rdm_kernel(h.orbs, h.pFrom, N, N)

    # TODO: optimize for sparse h
    H = np.zeros(K[0][0].shape)
    for i in range(len(K[0])):
        for j in range(len(K)):
            H += h.data[i, j] * K[j][i]

    return FermiOp(h.orbs, N, N, H)
