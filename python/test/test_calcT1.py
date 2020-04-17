import numpy as np
from scipy.special import binom
import fermifab
import unittest


class TestCalcT1(unittest.TestCase):

    def _err_calcT1(self, psi):
        g1 = fermifab.rdm(psi, 1)
        g2 = fermifab.rdm(psi, 2)

        # calculate T1 operator
        T = fermifab.calcT1(g2, g1)

        # T1 operator based solely on two-body RDM
        Td = fermifab.calcT1_(g2, psi.N)

        err = fermifab.norm(T - Td)

        # TODO: implement test for T1N

        # generalized T1 operator for p = 2 should be gamma^2(psi) + Q
        # TODO: implement test for T1 using p = 2

        return err

    def test_calcT1(self):
        psi = fermifab.FermiState(6, 4)
        x = fermifab.crand(int(binom(6, 4)))
        psi.data = x/np.linalg.norm(x)

        err = self._err_calcT1(psi)
        self.assertAlmostEqual(err, 0)


if __name__ == '__main__':
    unittest.main()