import numpy as np
from scipy.special import binom
import fermifab
import unittest


class TestNorbs(unittest.TestCase):

    def _norbs_err(self, orbs, N, real_valued=True):
        """norbs - 'natural orbitals': eigenstates of the 1-body RDM"""

        if real_valued:
            psi = fermifab.FermiState(orbs, N, data=np.random.rand(int(binom(orbs, N))))
        else:
            psi = fermifab.FermiState(orbs, N, data=fermifab.crand(int(binom(orbs, N))))

        # compute the one-body reduced density matrix
        rho = fermifab.rdm(psi, 1)

        # diagonalize rho and construct base change matrix operating on N-fold tensor product
        _, U = fermifab.eig(rho)
        UN = fermifab.tensor_op(U, N)
        # apply inverse base change matrix to psi
        psi = UN.H @ psi

        G = fermifab.rdm(psi, 1)

        # G should be a diagonal matrix
        return np.linalg.norm(np.diag(np.diag(G.data)) - G.data)


    def test_norbs(self):
        self.assertAlmostEqual(self._norbs_err(6, 4, True), 0)
        self.assertAlmostEqual(self._norbs_err(6, 4, False), 0)


if __name__ == '__main__':
    unittest.main()

