/// \file mextools.h
/// \brief Matlab 'mex' utility functions.
//
//	Copyright (c) 2008-2015, Christian B. Mendl
//	All rights reserved.
//	http://christian.mendl.net
//
//	This program is free software; you can redistribute it and/or
//	modify it under the terms of the Simplified BSD License
//	http://www.opensource.org/licenses/bsd-license.php
//
//	Reference:
//		Christian B. Mendl
//		The FermiFab toolbox for fermionic many-particle quantum systems
//		Comput. Phys. Commun. 182, 1327-1337 (2011)
//		preprint http://arxiv.org/abs/1103.0872
//________________________________________________________________________________________________________________________
//

#ifndef MEXTOOLS_H
#define MEXTOOLS_H

#include "util.h"
#include <mex.h>
#include <search.h>

// mwIndex unknown to older versions of Matlab
#if defined DEFINE_INDEX
typedef int mwIndex;
#endif


// create empty m x n sparse array
mxArray *CreateEmptySparse(const int m, const int n);


//________________________________________________________________________________________________________________________
///
/// \brief Matlab sparse array entry
///
typedef struct
{
	int i;			//!< zero-based index
	int x;			//!< data
}
sparse_array_entry_t;


// copy array holding real entries to Matlab sparse array
mxArray *ToSparseArray(const int n, sparse_array_entry_t *entries, int num);


//________________________________________________________________________________________________________________________
///
/// \brief Matlab sparse matrix entry
///
typedef struct
{
	int i, j;		//!< zero-based indices
	int x;			//!< data
}
sparse_matrix_entry_t;


//________________________________________________________________________________________________________________________
///
/// \brief Return true if indices are different
///
static inline bool DiffIndex(const sparse_matrix_entry_t *x, const sparse_matrix_entry_t *y)
{
	return (x->i != y->i) || (x->j != y->j);
}


// copy array holding real matrix entries to Matlab sparse matrix
mxArray *ToSparseMatrix(const int m, const int n, sparse_matrix_entry_t *entries, int num);



#endif
