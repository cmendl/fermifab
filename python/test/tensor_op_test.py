from scipy.special import binom
import numpy as np
import fermifab

def tensor_op_test(orbs, p, N):
    """TENSOR_OP_TEST - Operator tensor product test"""
    err = 0

    # TODO: allow complex data
    data = np.random.rand(int(binom(orbs, N)))
    psi = fermifab.FermiState(orbs, N, data = data)

    # compute one-body reduced density matrix
    rho = fermifab.rdm(psi, 1)
    U, _ = np.linalg.qr(rho.data)
    
    U = fermifab.FermiOp(orbs, 1, 1, U.real)

    Up = fermifab.tensor_op(U, p)
    UN = fermifab.tensor_op(U, N)

    err += np.linalg.norm((UN.T*UN).data - np.eye(*UN.shape))
    err += np.linalg.norm((Up.T*fermifab.rdm(UN*psi,p)*Up- fermifab.rdm(psi,p)).data)

    return err

if __name__ == '__main__':
    orbs = 6
    p = 2
    N = 4
    err = tensor_op_test(orbs, p, N)
    np.testing.assert_almost_equal(err, 0, 10)
