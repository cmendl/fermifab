from scipy.special import binom
import numpy as np
import fermifab
import unittest


class TestFermiState(unittest.TestCase):

    def test_arithmetic_operators(self):
        orbs = 6
        N = 4
        data1 = fermifab.crand(int(binom(orbs, N)))
        data2 = fermifab.crand(int(binom(orbs, N)))
        psi = fermifab.FermiState(orbs, N, data = data1)
        phi = fermifab.FermiState(orbs, N, data = data2)
        a, b = np.random.rand(2)

        err = np.linalg.norm((psi + phi).data - (psi.data+phi.data))
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((psi - phi).data - (psi.data-phi.data))
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((a*psi*b).data - a*psi.data*b)
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((psi/b).data - psi.data/b)
        self.assertAlmostEqual(err, 0)

        err = fermifab.norm(psi) - np.linalg.norm(psi.data)
        self.assertAlmostEqual(err, 0)


if __name__ == '__main__':
    unittest.main()