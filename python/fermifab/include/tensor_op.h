/// \file tensor_op.h
/// \brief Calculate the tensor product (A otimes A ... otimes A): wedge^N H -> wedge^N H for an operator A: H -> H.
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


double Det(const int n, double *A);


int TensorOp(const int orbs, const int N, const double *A, sparse_array_t *AN);
