/// \file gen_spintrace_sing.c
/// \brief Generate single-body spatial Hamiltonian mapping for 1 particle in 2*spatialorbs orbitals including spin,
/// i.e. the orbitals are (1 up, 1 down, 2 up, 2 down,...) with spatial orbitals 1, 2,...; Matlab 'mex' function.
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

#include "mextools.h"


static mxArray *GenerateSpintraceSingle(int spatialorbs);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_spintrace_sing(spatialorbs)
///
/// Generate single-body spatial Hamiltonian mapping; Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs > 1)  { mexErrMsgTxt("1 output argument expected."); }
	if (nrhs != 1) { mexErrMsgTxt("1 input arguments required."); }

	int spatialorbs = (int)mxGetScalar(prhs[0]);
	if (spatialorbs < 1) { mexErrMsgTxt("Number of spatial orbitals 'spatialorbs' must be positive."); }

	// verbose output
	// mexPrintf("Calculating spin-trace mapping for a single particle in %i orbitals (including spin)...", 2*spatialorbs);
	plhs[0] = GenerateSpintraceSingle(spatialorbs);
	// mexPrintf(" Done.\n");
}


//________________________________________________________________________________________________________________________
///
/// \brief Generate single-body spatial Hamiltonian mapping for 1 particle in 2*spatialorbs orbitals including spin,
/// i.e. the orbitals are (1 up, 1 down, 2 up, 2 down, ...) with spatial orbitals 1, 2, ...
///
static mxArray *GenerateSpintraceSingle(int spatialorbs)
{
	// allocate cell matrix array
	mxArray *ret = mxCreateCellMatrix(spatialorbs, spatialorbs);
	if (!ret) { mexErrMsgTxt("Creating cell matrix failed."); }

	sparse_matrix_entry_t matrix[2];
	matrix[0].x = 1;
	matrix[1].x = 1;

	mwIndex n[2];
	for (n[0] = 0; n[0] < (mwIndex)spatialorbs; n[0]++)
	{
		for (n[1] = 0; n[1] < (mwIndex)spatialorbs; n[1]++)
		{
			// both spins up
			matrix[0].j = (int)(2*n[0]);
			matrix[0].i = (int)(2*n[1]);

			// both spins down
			matrix[1].j = (int)(2*n[0]+1);
			matrix[1].i = (int)(2*n[1]+1);

			mxSetCell(ret, mxCalcSingleSubscript(ret, 2, n), ToSparseMatrix(2*spatialorbs, 2*spatialorbs, matrix, 2));
		}
	}

	return ret;
}
