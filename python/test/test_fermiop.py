from scipy.special import binom
import numpy as np
import fermifab
import unittest

class TestFermiOp(unittest.TestCase):
    def test_arithmetic_operators(self):
        orbs = 6
        N = 4
        data1 = fermifab.crand(int(binom(orbs, N)), int(binom(orbs, N)))
        data2 = fermifab.crand(int(binom(orbs, N)), int(binom(orbs, N)))
        S = fermifab.FermiOp(orbs, N, N, data = data1)
        R = fermifab.FermiOp(orbs, N, N, data = data2)
        a, b = fermifab.crand(2)

        err = np.linalg.norm((S + R).data - (S.data + R.data))
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((S - R).data - (S.data - R.data))
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((a*S*b).data - a*S.data*b)
        self.assertAlmostEqual(err, 0)

        err = np.linalg.norm((S/a).data - S.data/a)
        self.assertAlmostEqual(err, 0)

    def test_density_matrix(self):
        orbs = 6
        N = 4
        data = fermifab.crand(int(binom(orbs, N)))
        psi = fermifab.FermiState(orbs, N, data = data)
        psi = psi/fermifab.norm(psi)
        rho = psi @ psi.T

        
        self.assertAlmostEqual((psi.T @ psi),1)

        err = np.linalg.norm((rho @ psi).data - psi.data)
        self.assertAlmostEqual(err, 0)

        err = fermifab.norm(rho @ rho - rho)
        self.assertAlmostEqual(err, 0)
        

if __name__ == '__main__':
    unittest.main()