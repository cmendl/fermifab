import fermifab
import numpy as np
from scipy.special import binom

__all__ = ['calcQ', 'calcT1']

def calcQ(g2, g1):
    """
    Calculate the Q operator of the Q-representability condition
    from the one- and two-body RDMs g1 and g2, respectively

    If g1 is numeric, it is interpreted as the number of particles N
    and g1 is computed from g2. 

    Reference: "The reduced density matrix method for electronic structure
    calculations and the role of three-index representability conditions"
    Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721

    Copyright (c) 2008-2015, Christian B. Mendl
    All rights reserved.
    """

    if isinstance(g1, (float, complex, int)):
        N = g1
        g1 = fermifab.rdm(g2, 1)/(N-1)
    else:
        assert type(g1) == fermifab.FermiOp

    id = g2.copy()
    id.data = np.eye(*g2.shape)

    Q = g2 - fermifab.p2N(g1, 2) + id
    return Q

def calcT1(g2, g1):
    """
    Calculate the T1 operator of the T1-representability condition
    from the one- and two-body RDMs g1 and g2, respectively

    If g1 is numeric, it is interpreted as the number of particles N
    and g1 is computed from g2. 

    Reference: "The reduced density matrix method for electronic structure
    calculations and the role of three-index representability conditions"
    Journal of Chemical Physics 120, 2095 (2004); doi:10.1063/1.1636721

    Copyright (c) 2008-2015, Christian B. Mendl
    All rights reserved.
    """
    if isinstance(g1, (float, complex, int)):
        N = g1
        g1 = fermifab.rdm(g2, 1)/(N-1)
    else:
        assert type(g1) == fermifab.FermiOp

    id = fermifab.FermiOp(g2.orbs, 3, 3)
    id.data =  np.eye(int(binom(g2.orbs, 3)))

    T = fermifab.p2N(g2, 3) - fermifab.p2N(g1, 3) + id
    return T