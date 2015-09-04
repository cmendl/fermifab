/// \file mextools.c
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

#include <assert.h>
#include "mextools.h"


//________________________________________________________________________________________________________________________
///
/// \brief Create an empty m x n sparse array
///
mxArray *CreateEmptySparse(const int m, const int n)
{
	int j;

	// allocate sparse matrix
	mxArray *ret = mxCreateSparse(m, n, 0, mxREAL);

	// fill 'jc'
	mwIndex *jc = mxGetJc(ret);
	for (j = 0; j <= n; j++)
	{
		jc[j] = 0;
	}

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
static int CompareSparseArrayEntry(const void *x, const void *y)
{
	const sparse_array_entry_t *a = x;
	const sparse_array_entry_t *b = y;

	if (a->i < b->i) {
		return -1;
	}

	if (b->i < a->i) {
		return 1;
	}

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy array holding real entries to Matlab sparse array;
/// the values of multiply indexed entries are summed up
///
mxArray *ToSparseArray(const int n, sparse_array_entry_t *entries, int num)
{
	int j, k;

	// need entries ordered
	qsort(entries, num, sizeof(sparse_array_entry_t), CompareSparseArrayEntry);

	if (num > 0)
	{
		// sum up multiple entries

		j = 0;
		for (k = 1; k < num; k++)
		{
			if (entries[j].i != entries[k].i)
			{
				j++;
				if (j != k) {
					entries[j] = entries[k];
				}
			}
			else
			{
				entries[j].x += entries[k].x;
			}
		}
		num = j + 1;

		// remove zeros

		j = 0;
		for (k = 1; k < num; k++)
		{
			if (entries[k].x != 0)
			{
				j++;
				if (j != k) {
					entries[j] = entries[k];
				}
			}
		}
		num = j + 1;

		if (entries[0].x == 0)
		{
			// remove first entry
			entries++;
			num--;
		}
	}

	// integrity check: no multiple entries any more
	for (k = 1; k < num; k++)
	{
		assert(CompareSparseArrayEntry(&entries[k-1], &entries[k]) < 0);
	}

	// allocate sparse matrix
	mxArray *ret = mxCreateSparse(n, 1, num, mxREAL);

	// copy nonzeros
	double  *sr = mxGetPr(ret);
	mwIndex *ir = mxGetIr(ret);
	for (k = 0; k < num; k++)
	{
		sr[k] = entries[k].x;
		ir[k] = entries[k].i;
	}

	// fill 'jc'
	mwIndex *jc = mxGetJc(ret);
	jc[0] = 0;
	jc[1] = num;

	return ret;
}


//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
static int CompareSparseMatrixEntry(const void *x, const void *y)
{
	const sparse_matrix_entry_t *a = x;
	const sparse_matrix_entry_t *b = y;

	// sort column first

	if (a->j < b->j)
		return -1;

	if (b->j < a->j)
		return 1;

	// a->j == b->j

	if (a->i < b->i)
		return -1;

	if (b->i < a->i)
		return 1;

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Copy array holding real matrix entries to Matlab sparse matrix;
/// the values of multiply indexed entries are summed up
/// 
mxArray *ToSparseMatrix(const int m, const int n, sparse_matrix_entry_t *entries, int num)
{
	int j, k;

	// need entries ordered
	qsort(entries, num, sizeof(sparse_matrix_entry_t), CompareSparseMatrixEntry);

	if (num > 0)
	{
		// sum up multiple entries

		j = 0;
		for (k = 1; k < num; k++)
		{
			if (DiffIndex(&entries[j], &entries[k]))
			{
				j++;
				if (j != k)
				{
					entries[j] = entries[k];
				}
			}
			else
			{
				entries[j].x += entries[k].x;
			}
		}
		num = j + 1;

		// remove zeros

		j = 0;
		for (k = 1; k < num; k++)
		{
			if (entries[k].x != 0)
			{
				j++;
				if (j != k) {
					entries[j] = entries[k];
				}
			}
		}
		num = j + 1;

		if (entries[0].x == 0)
		{
			// remove first entry
			entries++;
			num--;
		}
	}

	// integrity check: no multiple entries any more
	for (k = 1; k < num; k++)
	{
		assert(CompareSparseMatrixEntry(&entries[k-1], &entries[k]) < 0);
	}

	// allocate sparse matrix
	mxArray *ret = mxCreateSparse(m, n, num, mxREAL);

	mwIndex *jc = mxGetJc(ret);
	if (num > 0)
	{
		// copy nonzeros
		double  *sr = mxGetPr(ret);
		mwIndex *ir = mxGetIr(ret);
		for (k = 0; k < num; k++)
		{
			sr[k] = entries[k].x;
			ir[k] = entries[k].i;
		}

		// fill 'jc'
		for (j = 0; j <= entries[0].j; j++) {
			jc[j] = 0;
		}
		for (k = 1; k < num; k++)
		{
			for (j = entries[k-1].j + 1; j <= entries[k].j; j++) {
				jc[j] = k;
			}
		}
		for (j = entries[num-1].j + 1; j <= n; j++) {
			jc[j] = num;
		}
	}
	else
	{
		// fill 'jc'
		for (j = 0; j <= n; j++) {
			jc[j] = 0;
		}
	}

	return ret;
}
