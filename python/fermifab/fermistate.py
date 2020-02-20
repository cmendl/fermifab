import numpy as np
from scipy.special import binom


__all__ = ['FermiState']


class FermiState(object):

    def __init__(self, orbs, N, data=None):
        self.orbs = orbs
        self.N = N
        if data != None:
            # TODO: support lists of 'orbs' and 'N'
            assert len(data) == int(binom(orbs, N))
            self.data = np.array(data)
        else:
            self.data = np.zeros(int(binom(orbs, N)), dtype=complex)
