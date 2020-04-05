from scipy.special import binom
import fermifab
import unittest


class TestP2N(unittest.TestCase):

    def _p2N_err(self, orbs, p, N):
        shape_data = int(binom(orbs, p)), int(binom(orbs, p))
        h = fermifab.FermiOp(orbs, p, p, fermifab.crand(*shape_data))

        # random wavefunction "psi"
        psi = fermifab.FermiState(orbs, N, data=fermifab.crand(int(binom(orbs, N))))

        return abs((psi.H @ fermifab.p2N(h, N) @ psi).data[0] - fermifab.trace(h @ fermifab.rdm(psi, p)))

    def test_p2N(self):
        self.assertAlmostEqual(self._p2N_err(6, 2, 4), 0)
        self.assertAlmostEqual(self._p2N_err(7, 1, 3), 0)


if __name__ == '__main__':
    unittest.main()
