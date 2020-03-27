import numpy as np
from scipy.special import binom
import fermifab
import fermifab.repr_conditions as RC
import unittest


class TestCalcQ(unittest.TestCase):
    def _err_calcQ(self, psi):
        g1 = fermifab.rdm(psi, 1)
        g2 = fermifab.rdm(psi, 2)
        
        # Calculate Q operator
        Q = RC.calcQ(g2, g1)
        
        ##
        # Q operator based solely on two-body RDM
        Q2 = RC.calcQ(g2, psi.N)
        err = fermifab.norm(Q - Q2)
        
        ## 
        # TODO: implement test for QN 

        ##
        # generalized Q operator for p = 1 should be id - gamma^1(psi)
        # TODO: implement test for QN using p = 1

        return err

    def test_calcQ(self):
        psi = fermifab.FermiState(6, 4, fermifab.crand(int(binom(6, 4))))

        err = self._err_calcQ(psi)
        self.assertAlmostEqual(err, 0)



if __name__ == '__main__':
    unittest.main()