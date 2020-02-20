/// \file sparse.h
/// \brief Sparse arrays.
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


//________________________________________________________________________________________________________________________
///
/// \brief Sparse array structure
///
typedef struct
{
	double *val;    //!< non-zero values
	int *ind;       //!< corresponding indices (matrix of dimension 'nnz x rank')
	int *dims;      //!< dimensions (vector of length 'rank')
	int nnz;        //!< number of non-zero entries
	int rank;       //!< array rank
}
sparse_array_t;


void DeleteSparseArray(sparse_array_t *a);
