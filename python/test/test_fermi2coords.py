from scipy.special import binom
import numpy as np
import fermifab
import unittest

class TestFermi2Coords(unittest.TestCase):
    def test_fermi2coords(self):
        X = fermifab.fermi2coords(10,4)
        X_matlab = np.loadtxt("test_fermi2coords.csv") 
        
        X_matlab += -1 # Matlab uses 1-based indexing

        err = np.linalg.norm(X - X_matlab)
        self.assertAlmostEqual(err, 0)
        
if __name__ == '__main__':
    unittest.main()