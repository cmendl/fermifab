/// \file calc_tensor_op.c
/// \brief Calculate the tensor product (A otimes A ... otimes A): wedge^N H -> wedge^N H for an operator A: H -> H; Matlab 'mex' function.
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

#include "fermi_map.h"
#include "util.h"
#include "mextools.h"


//________________________________________________________________________________________________________________________
///
/// \brief AN = calc_tensor_op(A,N)
///
/// Calculate the tensor product (A otimes A ... otimes A):
/// wedge^N H -> wedge^N H for an operator A: H -> H;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, k, l;

	if (nrhs != 2) { mexErrMsgTxt("2 input arguments required. "); }
	if (nlhs > 1)  { mexErrMsgTxt("1 output argument expected."); }

	const mxArray *A = prhs[0];
	int orbs = (int)mxGetM(A);
	int N    = (int)mxGetScalar(prhs[1]);

	// check input arguments
	if (mxIsSparse(A))     { mexErrMsgTxt("Input argument 'A' cannot be a sparse matrix."); }
	if (orbs != mxGetN(A)) { mexErrMsgTxt("Input argument 'A' must be a square matrix."); }
	if (N < 1)             { mexErrMsgTxt("Input argument 'N' must at least be 1."); }
	if (N > orbs)          { mexErrMsgTxt("Dimension of first input argument (i.e., 'orbs') must at least be 'N'."); }

	// create Fermi map
	fermi_map_t baseMap;
	{
		fermi_config_t config;
		config.orbs = &orbs;
		config.N    = &N;
		config.nc   = 1;
		FermiMap(&config, &baseMap);
	}

	fermi_coords_t *x = (fermi_coords_t *)mxMalloc(N*sizeof(fermi_coords_t));
	fermi_coords_t *y = (fermi_coords_t *)mxMalloc(N*sizeof(fermi_coords_t));

	if (mxIsComplex(A))
	{
		double *prA = mxGetPr(A);
		double *piA = mxGetPi(A);

		// temporary matrix for determinant calculation
		mxArray *T = mxCreateDoubleMatrix(N, N, mxCOMPLEX);

		// allocate return value matrix
		int nzmax = 16;
		plhs[0] = mxCreateSparse(baseMap.num, baseMap.num, nzmax, mxCOMPLEX);
		double *prAN = mxGetPr(plhs[0]);
		double *piAN = mxGetPi(plhs[0]);
		mwIndex *ir  = mxGetIr(plhs[0]);
		mwIndex *jc  = mxGetJc(plhs[0]);

		int nnz = 0;
		for (j = 0; j < baseMap.num; j++)
		{
			jc[j] = nnz;
			FermiDecode(baseMap.map[j], y, N);

			for (i = 0; i < baseMap.num; i++)
			{
				mxArray *det;

				// short-circuit zero determinant detection
				bool zero_det = false;

				double *prT = mxGetPr(T);
				double *piT = mxGetPi(T);

				FermiDecode(baseMap.map[i], x, N);

				// copy entries in 'A' indexed by x and y to 'T'
				for (l = 0; l < N; l++)
				{
					bool zero_col = true;
					for (k = 0; k < N; k++)
					{
						mwIndex indexA = (mwIndex)(x[k] + orbs*y[l]);
						*prT = prA[indexA];
						*piT = piA[indexA];
						if ((*prT) != 0 || (*piT) != 0) {
							zero_col = false;
						}
						prT++;
						piT++;
					}
					zero_det |= zero_col;
				}

				// short-circuit zero determinant detection
				if (zero_det) {
					continue;
				}

				if (mexCallMATLAB(1, &det, 1, &T, "det") != 0) {
					mexErrMsgTxt("Calling 'det' failed.");
				}

				if (*mxGetPr(det) == 0 && !mxIsComplex(det)) {
					continue;
				}

				// check whether we have to reallocate storage space
				if (nnz == nzmax)
				{
					nzmax += orbs;
					mxSetNzmax(plhs[0], nzmax); 
					mxSetPr(plhs[0],  (double *)mxRealloc(prAN, nzmax*sizeof(double)));
					mxSetPi(plhs[0],  (double *)mxRealloc(piAN, nzmax*sizeof(double)));
					mxSetIr(plhs[0], (mwIndex *)mxRealloc(ir,   nzmax*sizeof(mwIndex)));

					prAN = mxGetPr(plhs[0]);
					piAN = mxGetPi(plhs[0]);
					ir   = mxGetIr(plhs[0]);
				}

				prAN[nnz] = *mxGetPr(det);
				if (mxIsComplex(det)) {
					piAN[nnz] = *mxGetPi(det);
				}
				else {
					piAN[nnz] = 0;
				}
				ir[nnz] = i;

				nnz++;
			}
		}
		jc[baseMap.num] = nnz;

		mxDestroyArray(T);
	}
	else	// !mxIsComplex(A)
	{
		double *prA = mxGetPr(A);

		// temporary matrix for determinant calculation
		mxArray *T = mxCreateDoubleMatrix(N, N, mxREAL);

		// allocate return value matrix
		int nzmax = 16;
		plhs[0] = mxCreateSparse(baseMap.num, baseMap.num, nzmax, mxREAL);
		double *prAN = mxGetPr(plhs[0]);
		mwIndex *ir  = mxGetIr(plhs[0]);
		mwIndex *jc  = mxGetJc(plhs[0]);

		int nnz = 0;
		for (j = 0; j < baseMap.num; j++)
		{
			jc[j] = nnz;
			FermiDecode(baseMap.map[j], y, N);

			for (i = 0; i < baseMap.num; i++)
			{
				mxArray *det;

				// short-circuit zero determinant detection
				bool zero_det = false;

				double *prT = mxGetPr(T);

				FermiDecode(baseMap.map[i], x, N);

				// copy entries in 'A' indexed by x and y to 'T'
				for (l = 0; l < N; l++)
				{
					bool zero_col = true;
					for (k = 0; k < N; k++)
					{
						mwIndex indexA = (mwIndex)(x[k] + orbs*y[l]);
						*prT = prA[indexA];
						if ((*prT) != 0) {
							zero_col = false;
						}
						prT++;
					}
					zero_det |= zero_col;
				}

				// short-circuit zero determinant detection
				if (zero_det) {
					continue;
				}

				if (mexCallMATLAB(1, &det, 1, &T, "det") != 0) {
					mexErrMsgTxt("Calling 'det' failed.");
				}

				if (*mxGetPr(det) == 0) {
					continue;
				}

				// check whether we have to reallocate storage space
				if (nnz == nzmax)
				{
					nzmax += orbs;
					mxSetNzmax(plhs[0], nzmax); 
					mxSetPr(plhs[0],  (double *)mxRealloc(prAN, nzmax*sizeof(double)));
					mxSetIr(plhs[0], (mwIndex *)mxRealloc(ir,   nzmax*sizeof(mwIndex)));

					prAN = mxGetPr(plhs[0]);
					ir   = mxGetIr(plhs[0]);
				}

				prAN[nnz] = *mxGetPr(det);
				ir[nnz] = i;

				nnz++;
			}
		}
		jc[baseMap.num] = nnz;

		mxDestroyArray(T);
	}

	// clean up
	mxFree(y);
	mxFree(x);
	mxFree(baseMap.map);
}
