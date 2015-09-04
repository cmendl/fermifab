/// \file gen_annihil.c
/// \brief Generate matrix representation of an annihilation operator; Matlab 'mex' function.
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
#include "mextools.h"


static mxArray *GenerateAnnihilation(const int orbs, const int p, const int N);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_annihil(orbs,p,N)
///
/// Generate matrix representation of an annihilation operator;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3) { mexErrMsgTxt("3 input arguments required."); }
	if (nlhs > 1)  { mexErrMsgTxt("1 output argument expected."); }

	int orbs = (int)mxGetScalar(prhs[0]);
	int p    = (int)mxGetScalar(prhs[1]);
	int N    = (int)mxGetScalar(prhs[2]);

	if (p < 1)    { mexErrMsgTxt("'p' must at least be 1."); }
	if (p >= N)   { mexErrMsgTxt("'p' must be less than 'N'."); }
	if (N > orbs) { mexErrMsgTxt("Number of orbitals 'orbs' must at least be 'N'."); }

	// verbose output
	// mexPrintf("Creating annihilation operator (p = %i, N = %i, number of orbitals: %i)...", p, N, orbs);
	plhs[0] = GenerateAnnihilation(orbs, p, N);
	// mexPrintf(" Done.\n");
}


//________________________________________________________________________________________________________________________
///
/// \brief Generate matrix representation of an annihilation operator on a Hilbert space with dimension 'orbs',
/// i.e., given a normalized p-body wavefunction chi (coefficients of canonical Slater basis),
/// chi' K{i,j} is the coefficient i,j of a_chi operating on wedge^N H
///
static mxArray *GenerateAnnihilation(const int orbs, const int p, const int N)
{
	// Fermi map with 'p' particles
	fermi_map_t baseMapP;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbs;
		config.N    = (int *)&p;
		config.nc   = 1;
		FermiMap(&config, &baseMapP);
	}
	// Fermi map with 'N' particles
	fermi_map_t baseMapN;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbs;
		config.N    = (int *)&N;
		config.nc   = 1;
		FermiMap(&config, &baseMapN);
	}
	// Fermi map with 'N - p' particles
	const int Np = N - p;
	fermi_map_t baseMapNP;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbs;
		config.N    = (int *)&Np;
		config.nc   = 1;
		FermiMap(&config, &baseMapNP);
	}

	// allocate cell matrix array
	mxArray *ret = mxCreateCellMatrix(baseMapNP.num, baseMapN.num);
	if (!ret) { mexErrMsgTxt("Creating cell matrix failed (probably out of memory)."); }

	// initialize with zeros
	mwIndex n[2];
	for (n[0] = 0; n[0] < (mwIndex)baseMapNP.num; n[0]++)
	{
		for (n[1] = 0; n[1] < (mwIndex)baseMapN.num; n[1]++)
		{
			mxSetCell(ret, mxCalcSingleSubscript(ret, 2, n), CreateEmptySparse(baseMapP.num, 1));
		}
	}

	fermi_coords_t *x = (fermi_coords_t *)mxMalloc(N*sizeof(fermi_coords_t));

	sparse_array_entry_t entry;
	for (entry.i = 0; entry.i < baseMapP.num; entry.i++)
	{
		FermiDecode(baseMapP.map[entry.i], x, p);

		for (n[0] = 0; n[0] < (mwIndex)baseMapNP.num; n[0]++)
		{
			FermiDecode(baseMapNP.map[n[0]], x + p, Np);

			n[1] = Fermi2BaseSign(&baseMapN, x, N, &entry.x);
			if (n[1] == -1) {
				continue;
			}

			mxSetCell(ret, mxCalcSingleSubscript(ret, 2, n), ToSparseArray(baseMapP.num, &entry, 1));
		}
	}

	// clean up
	mxFree(x);
	mxFree(baseMapNP.map);
	mxFree(baseMapN.map);
	mxFree(baseMapP.map);

	return ret;
}
