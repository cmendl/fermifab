/// \file sparse.c
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

#include "sparse.h"
#include "util.h"
#include <stdlib.h>
#include <memory.h>


void DeleteSparseArray(sparse_array_t *a)
{
	if (a->dims != NULL) {  free(a->dims); }
	if (a->ind  != NULL) {  free(a->ind);  }
	if (a->val  != NULL) {  free(a->val);  }

	a->nnz  = 0;
	a->rank = 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Convert tensor index to data offset (row major ordering)
///
static inline int IndexToOffset(const int rank, const int *dims, const int *ind)
{
	int i;
	int offset = 0;
	int stride = 1;
	for (i = rank-1; i >= 0; i--)
	{
		offset += stride*ind[i];
		stride *= dims[i];
	}

	return offset;
}


//________________________________________________________________________________________________________________________
///
/// \brief Convert a sparse to a dense array, assuming that memory for target array has been allocated
///
void SparseToDense(const sparse_array_t *a, double *dns)
{
	memset(dns, 0, IntegerProduct(a->dims, a->rank) * sizeof(double));

	int i;
	for (i = 0; i < a->nnz; i++)
	{
		dns[IndexToOffset(a->rank, a->dims, &a->ind[a->rank*i])] = a->val[i];
	}
}
