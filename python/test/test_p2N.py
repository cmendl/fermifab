import numpy as np
from scipy.special import binom
import fermifab
import unittest


class TestP2N(unittest.TestCase):
    def _p2N_err(self, orbs, p, N):
        shape_data = int(binom(orbs, p)), int(binom(orbs, p))
        h = fermifab.FermiOp(orbs, p, p, fermifab.crand(*shape_data))

        # random wavefunction "psi"
        psi = fermifab.FermiState(orbs, N, data=fermifab.crand(int(binom(orbs, N))))
    
        err = abs(psi.H@fermifab.p2N(h, N)@psi - fermifab.trace(h@fermifab.rdm(psi, p)))
        return err
    
    def test_p2N(self):
        err = self._p2N_err(6, 2, 4)
        self.assertAlmostEqual(err, 0)



if __name__ == '__main__':
    unittest.main()