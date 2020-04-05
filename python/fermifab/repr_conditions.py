import numpy as np
from scipy.special import binom
from .fermiop import FermiOp
from .rdm import rdm
from .p2N import p2N

__all__ = ['calcQ', 'calcQ_', 'calcT1', 'calcT1_']


def calcQ(g2, g1):
    """
    Calculate the Q operator of the Q-representability condition
    from the one- and two-body RDMs `g1` and `g2`, respectively.

    Reference:
        Z. Zhao, B. J. Braams, M. Fukuda, M. L. Overton, J. K. Percus
        The reduced density matrix method for electronic structure calculations and the role of three-index representability conditions
        J. Chem. Phys. 120, 2095 (2004); doi:10.1063/1.1636721
    """
    idop = g2.copy()
    idop.data = np.eye(*g2.shape)

    return g2 - p2N(g1, 2) + idop


def calcQ_(g2, N):
    """
    Calculate the Q operator of the Q-representability condition
    from the two-body RDM `g2`, with `N` the overall particle number.

    Reference:
        Z. Zhao, B. J. Braams, M. Fukuda, M. L. Overton, J. K. Percus
        The reduced density matrix method for electronic structure calculations and the role of three-index representability conditions
        J. Chem. Phys. 120, 2095 (2004); doi:10.1063/1.1636721
    """
    g1 = rdm(g2, 1)/(N - 1)
    return calcQ(g2, g1)


def calcT1(g2, g1):
    """
    Calculate the T1 operator of the T1-representability condition
    from the one- and two-body RDMs `g1` and `g2`, respectively.

    Reference:
        Z. Zhao, B. J. Braams, M. Fukuda, M. L. Overton, J. K. Percus
        The reduced density matrix method for electronic structure calculations and the role of three-index representability conditions
        J. Chem. Phys. 120, 2095 (2004); doi:10.1063/1.1636721
    """
    idop = FermiOp(g2.orbs, 3, 3)
    idop.data =  np.eye(int(binom(g2.orbs, 3)))

    return p2N(g2, 3) - p2N(g1, 3) + idop


def calcT1_(g2, N):
    """
    Calculate the T1 operator of the T1-representability condition
    from the two-body RDM `g2`, with `N` the overall particle number.

    Reference:
        Z. Zhao, B. J. Braams, M. Fukuda, M. L. Overton, J. K. Percus
        The reduced density matrix method for electronic structure calculations and the role of three-index representability conditions
        J. Chem. Phys. 120, 2095 (2004); doi:10.1063/1.1636721
    """
    g1 = rdm(g2, 1)/(N - 1)
    return calcT1(g2, g1)
