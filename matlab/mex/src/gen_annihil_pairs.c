/// \file gen_annihil_pairs.c
/// \brief Generate matrix representation of an annihilation operator, restricted to pair structure wavefunctions; Matlab 'mex' function.
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


static mxArray *GenerateAnnihilationPairs(const int orbshalf, const int phalf, const int Nhalf);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_annihil_pairs(orbshalf,phalf,Nhalf)
///
/// Generate matrix representation of an annihilation operator, but restricted to pair structure wavefunctions;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 3) { mexErrMsgTxt("3 input arguments required."); }
	if (nlhs > 1)  { mexErrMsgTxt("Wrong number of output arguments."); }

	int orbshalf = (int)mxGetScalar(prhs[0]);
	int phalf    = (int)mxGetScalar(prhs[1]);
	int Nhalf    = (int)mxGetScalar(prhs[2]);

	if (phalf < 1)        { mexErrMsgTxt("'phalf' must at least be 1."); }
	if (phalf >= Nhalf)   { mexErrMsgTxt("'phalf' must be less than 'Nhalf'."); }
	if (Nhalf > orbshalf) { mexErrMsgTxt("Number of orbitals 'orbshalf' must at least be 'Nhalf'."); }

	// verbose output
	// mexPrintf("Creating pair annihilation operator (p = %i, N = %i, number of orbitals: %i)...", 2*phalf, 2*Nhalf, 2*orbshalf);
	plhs[0] = GenerateAnnihilationPairs(orbshalf, phalf, Nhalf);
	// mexPrintf(" Done.\n");
}


//________________________________________________________________________________________________________________________
///
/// \brief Generate matrix representation of an annihilation operator on a Hilbert space,
/// similar to 'GenerateAnnihilation', but restricted to pair structure wavefunctions;
/// effectively just disregarding permutation signs
///
static mxArray *GenerateAnnihilationPairs(const int orbshalf, const int phalf, const int Nhalf)
{
	// Fermi map with 'phalf' particles
	fermi_map_t baseMapP;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbshalf;
		config.N    = (int *)&phalf;
		config.nc   = 1;
		FermiMap(&config, &baseMapP);
	}
	// Fermi map with 'Nhalf' particles
	fermi_map_t baseMapN;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbshalf;
		config.N    = (int *)&Nhalf;
		config.nc   = 1;
		FermiMap(&config, &baseMapN);
	}
	// Fermi map with 'Nhalf - phalf' particles
	const int Nphalf = Nhalf - phalf;
	fermi_map_t baseMapNP;
	{
		fermi_config_t config;
		config.orbs = (int *)&orbshalf;
		config.N    = (int *)&Nphalf;
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

	fermi_coords_t *x = (fermi_coords_t *)mxMalloc(Nhalf*sizeof(fermi_coords_t));

	sparse_array_entry_t entry;
	entry.x = 1;	// set "sign" to 1
	for (entry.i = 0; entry.i < baseMapP.num; entry.i++)
	{
		FermiDecode(baseMapP.map[entry.i], x, phalf);

		for (n[0] = 0; n[0] < (mwIndex)baseMapNP.num; n[0]++)
		{
			FermiDecode(baseMapNP.map[n[0]], x + phalf, Nphalf);

			n[1] = Fermi2Base(&baseMapN, x, Nhalf);
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
