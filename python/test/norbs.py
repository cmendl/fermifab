from scipy.special import binom
import numpy as np
import fermifab

def norbs(orbs, N):
    """norbs - 'natural orbitals': eigenstates of the 1-body RDM"""
    # TODO: allow complex data
    data = np.random.rand(int(binom(orbs, N)))
    psi = fermifab.FermiState(orbs, N, data = data)

    # compute one-body reduced density matrix
    rho = fermifab.rdm(psi, 1)
    _, U = np.linalg.eig(rho.data)
    
    U = U.real 

    state = fermifab.FermiOp(orbs, 1, 1, data=U)
    stateN = fermifab.tensor_op(state, N)
    psi = stateN.T * psi

    G = fermifab.rdm(psi,1).data

    err = np.linalg.norm(np.diag(np.diag(G)) - G)
    return psi, err

if __name__ == '__main__':
    orbs = 6
    N = 4
    _, err = norbs(orbs, N)
    np.testing.assert_almost_equal(err, 0, 10)
