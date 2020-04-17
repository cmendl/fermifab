from scipy.special import binom
import numpy as np
import fermifab
import unittest


class TestTensorOp(unittest.TestCase):

    def _tensor_op_err(self, orbs, p, N):
        err = 0

        # random wavefunction "psi"
        psi = fermifab.FermiState(orbs, N, data=fermifab.crand(int(binom(orbs, N))))

        # compute one-body reduced density matrix
        rho = fermifab.rdm(psi, 1)
        U, _ = np.linalg.qr(rho.data)

        U = fermifab.FermiOp(orbs, 1, 1, U)

        Up = fermifab.tensor_op(U, p)
        UN = fermifab.tensor_op(U, N)

        err += np.linalg.norm((UN.H @ UN).data - np.eye(*UN.shape))
        err += np.linalg.norm((Up.H @ fermifab.rdm(UN @ psi, p) @ Up - fermifab.rdm(psi, p)).data)
        return err

    def test_tensor_op(self):
        orbs = 6
        p = 2
        N = 4
        err = self._tensor_op_err(orbs, p, N)
        self.assertAlmostEqual(err, 0)


if __name__ == '__main__':
    unittest.main()
