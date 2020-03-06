/// \file tensor_op.c
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

#include "tensor_op.h"
#include "fermi_map.h"
#include "util.h"
#include <lapacke.h>
#include <stdbool.h>
#include <math.h>
#include <malloc.h>
#include <memory.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Compute the determinant of a real matrix by LU decomposition.
///
/// Input matrix `A` will be overwritten.
///
double Det(const int n, double *A)
{
	lapack_int *ipiv = (lapack_int *)malloc(n * sizeof(lapack_int));
	lapack_int info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, A, n, ipiv);
	if (info < 0)
	{
		free(ipiv);
		return NAN;
	}

	int sign = 1;
	int i;
	for (i = 0; i < n; i++)
	{
		// 'ipiv' uses 1-based indexing!
		if (ipiv[i] != i + 1) {
			sign = -sign;
		}
	}

	double logdet = 0;
	for (i = 0; i < n; i++)
	{
		if (A[i*(n + 1)] < 0) {
			sign = -sign;
		}
		logdet += log(fabs(A[i*(n + 1)]));
	}

	// clean up
	free(ipiv);

	return sign * exp(logdet);
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the tensor product (A otimes A ... otimes A): wedge^N H -> wedge^N H for an operator A: H -> H
///
int TensorOp(const int orbs, const int N, const double *A, sparse_array_t *AN)
{
	int i;
	int status;

	// create Fermi map
	fermi_map_t baseMap;
	{
		fermi_config_t config;
		config.orbs = &orbs;
		config.N    = &N;
		config.nc   = 1;
		status = FermiMap(&config, &baseMap);
		if (status < 0) { return status; }
	}

	fermi_coords_t *x = (fermi_coords_t *)malloc(N*sizeof(fermi_coords_t));
	fermi_coords_t *y = (fermi_coords_t *)malloc(N*sizeof(fermi_coords_t));

	// temporary matrix for determinant calculation
	double *T = (double *)malloc(N*N * sizeof(double));
	if (T == NULL) { return -1; }

	// setup to-be returned sparse matrix
	AN->rank = 2;
	AN->dims = (int *)malloc(AN->rank * sizeof(int));
	AN->dims[0] = baseMap.num;
	AN->dims[1] = baseMap.num;
	AN->nnz = 0;
	int nzmax = 16;
	AN->val = (double *)malloc(nzmax * sizeof(double));     if (AN->val == NULL) { return -1; }
	AN->ind = (int *)malloc(AN->rank*nzmax * sizeof(int));  if (AN->ind == NULL) { return -1; }

	for (i = 0; i < baseMap.num; i++)
	{
		FermiDecode(baseMap.map[i], x, N);

		int j;
		for (j = 0; j < baseMap.num; j++)
		{
			FermiDecode(baseMap.map[j], y, N);

			// short-circuit zero determinant detection
			bool zero_det = false;

			// copy entries in 'A' indexed by x and y to 'T'
			int k;
			for (k = 0; k < N; k++)
			{
				bool zero_row = true;
				int l;
				for (l = 0; l < N; l++)
				{
					T[N*k + l] = A[orbs*x[k] + y[l]];
					if (T[N*k + l] != 0) {
						zero_row = false;
					}
				}
				zero_det |= zero_row;
			}

			// short-circuit zero determinant detection
			if (zero_det) {
				continue;
			}

			double d = Det(N, T);
			if (d == 0) {
				continue;
			}

			// check whether we have to reallocate storage space
			if (AN->nnz == nzmax)
			{
				nzmax += orbs;
				AN->val = (double *)realloc(AN->val, nzmax * sizeof(double));     if (AN->val == NULL) { return -1; }
				AN->ind = (int *)realloc(AN->ind, AN->rank*nzmax * sizeof(int));  if (AN->ind == NULL) { return -1; }
			}

			// set entry
			AN->val[AN->nnz] = d;
			AN->ind[2*AN->nnz  ] = i;
			AN->ind[2*AN->nnz+1] = j;

			AN->nnz++;
		}
	}

	// clean up
	free(T);
	free(y);
	free(x);
	free(baseMap.map);

	return 0;
}
