import numpy as np
import fermifab
import unittest


class TestFermi2Coords(unittest.TestCase):

    def test_fermi2coords(self):

        X = fermifab.kernel.fermi2coords([5, 4], [2, 1])

        X_ref = np.array([
            [0, 1, 5], [0, 2, 5], [1, 2, 5], [0, 3, 5], [1, 3, 5],
            [2, 3, 5], [0, 4, 5], [1, 4, 5], [2, 4, 5], [3, 4, 5],
            [0, 1, 6], [0, 2, 6], [1, 2, 6], [0, 3, 6], [1, 3, 6],
            [2, 3, 6], [0, 4, 6], [1, 4, 6], [2, 4, 6], [3, 4, 6],
            [0, 1, 7], [0, 2, 7], [1, 2, 7], [0, 3, 7], [1, 3, 7],
            [2, 3, 7], [0, 4, 7], [1, 4, 7], [2, 4, 7], [3, 4, 7],
            [0, 1, 8], [0, 2, 8], [1, 2, 8], [0, 3, 8], [1, 3, 8],
            [2, 3, 8], [0, 4, 8], [1, 4, 8], [2, 4, 8], [3, 4, 8]])

        err = np.linalg.norm(X - X_ref)
        self.assertEqual(err, 0)


if __name__ == '__main__':
    unittest.main()
