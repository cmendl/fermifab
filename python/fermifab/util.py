import numpy as np
import fermifab

# TODO: Export kron when lists of orbs and N are implemented
__all__ = ['crand','norm','trace','matrix_power','comprise_config', 'eig']

def crand(*args):
    return 0.5 - np.random.rand(*args) + 1j*(0.5 - np.random.rand(*args))

def norm(x):
    if type(x) in [fermifab.FermiState,fermifab.FermiOp]:
        return float(np.linalg.norm(x.data))
    else:
        return TypeError("Argument must be FermiState or FermiOp")

def kron(x, y):
    """ Kronecker tensor product of two Fermi operators """
    if type(x) == type(y) == fermifab.FermiOp:
        return fermifab.FermiOp([x.orbs, y.orbs], 
                                [x.pFrom, y.pFrom],
                                [x.pTo,   y.pTo],
                                np.kron(x.data, y.data))

    elif type(x) == type(y) == fermifab.FermiState:
        return fermifab.FermiState([x.orbs, y.orbs],
                                   [x.N,    y.N],
                                   np.kron(x.data, y.data))
    else:
        raise TypeError("x and y must be both FermiOp or FermiState.")
    

def matrix_power(x, n):
    """ Matrix power of a Fermi operator """
    assert x.pFrom == x.pTo
    data = np.linalg.matrix_power(x.data, n)
    return fermifab.FermiOp(x.orbs, x.pFrom, x.pTo, data = data)

def eig(a):
    """Eigenvalues and eigenstates of Fermi operators"""
    assert type(a) == fermifab.FermiOp
    D, U = np.linalg.eig(a.data)
    return D, fermifab.FermiOp(a.orbs, a.pFrom, a.pTo, U)

def trace(a):
    """ Trace of a FermiOp"""
    assert a.pFrom == a.pTo
    return np.trace(a.data)

def comprise_config(orbs1, orbs2, N1, N2):
    """Comprise configurations"""

    assert (sum(orbs1) == sum(orbs2)) & (sum(N1) == sum(N2))
    assert np.all(orbs1 > 0) & np.all(orbs2>0)

    # cumulative number of orbitals
    cumorbs1 = np.cumsum(orbs1)
    cumorbs2 = np.cumsum(orbs2)
    #cumulative particle numbers
    n1 = np.cumsum(N1)
    n2 = np.cumsum(N2)

    # merge tables
    m1 = np.ones_like(cumorbs1,dtype = bool)
    m2 = np.ones_like(cumorbs2,dtype = bool)
    orbs = [orbs1[0]]
    N = [N1[0]]
    j1 = 0
    j2 = 0
    while j1 < len(cumorbs1)-1:
        if cumorbs1[j1] < cumorbs2[j2]:
            orbs[-1] += orbs1[j1+1]
            N[-1] += N1[j1+1]
            j1 += 1
        elif cumorbs2[j2] < cumorbs1[j1]:
            j2 += 1
        else:
            if n1[j1] == n2[j2]:
                # both mode and particle numbers match
                m1[j1] = False
                m2[j2] = False
                orbs.append(orbs1[j1+1])
                N.append(N1[j1+1])
            else:
                orbs[-1] += orbs1[j1+1]
                N[-1] += N1[j1+1]
            j1 += 1
            j2 += 1
    assert (sum(orbs) == cumorbs1[-1]) & (sum(N) == n1[-1])
    m1 = m1[:-1]
    m2 = m2[:-1]

    return m1, m2, orbs, N
