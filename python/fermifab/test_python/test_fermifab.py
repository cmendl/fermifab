from scipy.special import binom
import numpy as np
import fermifab

def test_fermistate():
    orbs = 6
    N = 4
    data = np.random.rand(int(binom(orbs, N)))
    psi = fermifab.FermiState(orbs, N, data = data)

    np.testing.assert_equal((psi/2+psi).data,psi.data/2+psi.data)
    np.testing.assert_equal((2*psi-psi).data,2*psi.data-psi.data)

def test_fermiop():
    orbs = 6
    N = 4
    data = np.random.rand(int(binom(orbs, N)))
    psi = fermifab.FermiState(orbs, N, data = data)
    psi = psi/fermifab.norm(psi)
    rho = psi @ psi.T

    np.testing.assert_almost_equal((psi.T @ psi),1)
    np.testing.assert_almost_equal((rho @ psi).data, psi.data)
    

if __name__ == "__main__":
    test_fermistate()
    test_fermiop()