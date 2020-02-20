/// \file generate_rdm.h
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices.
//
//  Copyright (c) 2008-2020, Christian B. Mendl
//  All rights reserved.
//  http://christian.mendl.net
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the Simplified BSD License
//  http://www.opensource.org/licenses/bsd-license.php
//
//  Reference:
//      Christian B. Mendl
//      The FermiFab toolbox for fermionic many-particle quantum systems
//      Comput. Phys. Commun. 182, 1327-1337 (2011)
//      preprint http://arxiv.org/abs/1103.0872
//________________________________________________________________________________________________________________________
//

#pragma once

#include "sparse.h"


int GenerateRDM(const int *orbs, const int *p1, const int *N1, const int *N2, const int nc, sparse_array_t *K);
