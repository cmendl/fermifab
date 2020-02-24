import numpy as np
from scipy.special import binom


__all__ = ['FermiOp']


class FermiOp(object):

    def __init__(self, orbs, pFrom, pTo, data=None):
        self.orbs = orbs
        self.pFrom = pFrom
        self.pTo = pTo
        H, W = (int(binom(orbs, pTo)), int(binom(orbs, pFrom)))
        if data is not None:
            assert data.shape == (H,W)
            self.data = np.array(data)
        else:
            self.data = np.zeros((H,W), dtype=complex)

