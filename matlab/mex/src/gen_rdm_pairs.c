/// \file gen_rdm_pairs.c
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices, restricted to pair structure wavefunctions; Matlab 'mex' function.
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
#include <assert.h>


static mxArray *GenerateRDMpairs(const int *orbshalf, const int *p1half, const int *N1half, const int *N2half, const int nc);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_rdm_pairs(orbshalf,phalf,Nhalf) \n
///        K = gen_rdm_pairs(orbshalf,phalf,N1half,N2half)
///
/// Create the matrix kernels K required for calculating p-body reduced density matrices,
/// restricted to pair structure wavefunctions;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;

	if (nrhs < 3 || nrhs > 4) { mexErrMsgTxt("3 or 4 input arguments required."); }
	if (nlhs > 1) { mexErrMsgTxt("1 output argument expected."); }

	// number of orbitals
	int nc;
	if (mxGetM(prhs[0]) == 1)
	{
		nc = (int)mxGetN(prhs[0]);
	}
	else
	{
		nc = (int)mxGetM(prhs[0]);
		if (mxGetN(prhs[0]) != 1) {
			mexErrMsgTxt("'orbshalf' input argument must be a vector of positive integers.");
		}
	}
	if (nc == 0) { mexErrMsgTxt("No 'orbshalf' input argument."); }
	int *orbshalf = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		orbshalf[i] = (int)(*(mxGetPr(prhs[0])+i));
		if (orbshalf[i] < 1) {
			mexErrMsgTxt("Number of orbitals must be positive.");
		}
	}

	// reduced particle configuration
	if ((mxGetM(prhs[1])!=nc || mxGetN(prhs[1])!=1) && (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=nc)) {
		mexErrMsgTxt("Number of orbitals and reduced particle fields must have the same length.");
	}
	int phalfTot = 0;	// half of total reduced particle number
	int *p1half = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p1half[i] = (int)(*(mxGetPr(prhs[1])+i));
		if (p1half[i] < 0) { mexErrMsgTxt("No negative particle numbers allowed."); }
		phalfTot += p1half[i];
	}
	if (phalfTot < 1) { mexErrMsgTxt("Reduced particle number must at least be 1."); }

	// particle configurations of full fermionic spaces
	if ((mxGetM(prhs[2])!=nc || mxGetN(prhs[2])!=1) && (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=nc)) {
		mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
	}
	int NhalfTot = 0;	// half of total particle number
	int *N1half = (int *)mxMalloc(nc*sizeof(int));
	int *N2half;
	for (i = 0; i < nc; i++)
	{
		N1half[i] = (int)(*(mxGetPr(prhs[2])+i));
		if (N1half[i] < 0 || N1half[i] > orbshalf[i]) {
			mexErrMsgTxt("Particle numbers must be between zero and the corresponding number of orbitals.");
		}
		NhalfTot += N1half[i];
	}
	if (nrhs == 4)
	{
		if ((mxGetM(prhs[3])!=nc || mxGetN(prhs[3])!=1) && (mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=nc)) {
			mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
		}

		int Ndiffhalf = 0;
		N2half = (int *)mxMalloc(nc*sizeof(int));
		for (i = 0; i < nc; i++)
		{
			N2half[i] = (int)(*(mxGetPr(prhs[3])+i));
			if (N2half[i] < 0 || N2half[i] > orbshalf[i]) {
				mexErrMsgTxt("Particle numbers must be between zero and the corresponding number of orbitals.");
			}
			Ndiffhalf += N1half[i] - N2half[i];
		}
		// check whether total particle numbers agree
		if (Ndiffhalf != 0) {
			mexErrMsgTxt("Total particle numbers in 'N1half' and 'N2half' must be the same.");
		}
	}
	else
	{
		N2half = N1half;
	}

	// consistency checks
	if (phalfTot > NhalfTot) { mexErrMsgTxt("Reduced particle number must be less or equal to the number of particles in the full system."); }
	for (i = 0; i < nc; i++)
	{
		if (p1half[i] > N2half[i]) {
			mexErrMsgTxt("'p1half' must be less or equal to the number of particles in the full system.");
		}
		if (N1half[i] - N2half[i] + p1half[i] < 0)	// p2half[i]
		{
			mexErrMsgTxt("No negative particle numbers allowed for p2 := N1 - N2 + p1.");
		}
	}

	// verbose output
	// mexPrintf("Creating reduced pair %i-body matrices from %i-fermion states (number of orbitals: %i)...", 2*phalfTot, 2*NhalfTot, 2*IntegerSum(orbshalf, nc));
	plhs[0] = GenerateRDMpairs(orbshalf, p1half, N1half, N2half, nc);
	// mexPrintf(" Done.\n");

	// clean up
	mxFree(N1half);
	if (N1half != N2half) {
		mxFree(N2half);
	}
	mxFree(p1half);
	mxFree(orbshalf);
}


//________________________________________________________________________________________________________________________
///
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices
/// similar to 'GenerateRdm', but restricted to pair structure wavefunctions,
/// i.e. |2i-1,2i> is treated as a virtual "bosonic" orbital
///
static mxArray *GenerateRDMpairs(const int *orbshalf, const int *p1half, const int *N1half, const int *N2half, const int nc)
{
	int i;

	const int N1halfTot = IntegerSum(N1half, nc);
	const int N2halfTot = IntegerSum(N2half, nc);
	const int p1halfTot = IntegerSum(p1half, nc);
	const int p2halfTot =  N1halfTot - N2halfTot + p1halfTot;

	// basic configuration setup
	fermi_config_t config;
	config.orbs = (int *)orbshalf;
	config.nc = nc;

	// fermionic particle spaces with configurations N1half and N2half
	fermi_map_t baseMapN1, baseMapN2;
	config.N = (int *)N1half; FermiMap(&config, &baseMapN1);
	config.N = (int *)N2half; FermiMap(&config, &baseMapN2);

	// configurations p1half and p2half
	fermi_map_t baseMapP1, baseMapP2;
	config.N = (int *)p1half; FermiMap(&config, &baseMapP1);
	// N2 - p1 == N1 - p2
	int *p2half = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p2half[i] = N1half[i] - N2half[i] + p1half[i];
		assert(p2half[i] >= 0);
	}
	config.N = p2half; FermiMap(&config, &baseMapP2);

	// allocate cell matrix array
	mxArray *ret = mxCreateCellMatrix(baseMapP1.num, baseMapP2.num);
	if (!ret) { mexErrMsgTxt("Creating cell matrix failed (probably out of memory)."); }

	fermi_coords_t *y = (fermi_coords_t *)mxMalloc(N1halfTot*sizeof(fermi_coords_t));

	mwIndex n[2];
	for (n[0] = 0; n[0] < (mwIndex)baseMapP1.num; n[0]++)
	{
		for (n[1] = 0; n[1] < (mwIndex)baseMapP2.num; n[1]++)
		{
			// fill first 'p2halfTot' coordinates of 'y' only
			FermiDecode(baseMapP2.map[n[1]], y, p2halfTot);

			// determine number of matrix entries for current p1 and p2 configuration
			int numEntries = 1;
			bitfield_t conf = baseMapP1.map[n[0]] | baseMapP2.map[n[1]];
			for (i = 0; i < nc; i++)
			{
				bitfield_t mask = ((bitfield_t)1 << orbshalf[i]) - 1;
				numEntries *= Binomial(orbshalf[i] - BitCount(mask & conf), N2half[i] - p1half[i]);
				conf <<= orbshalf[i];
			}

			sparse_matrix_entry_t *matrix = (sparse_matrix_entry_t *)mxMalloc(numEntries*sizeof(sparse_matrix_entry_t));

			sparse_matrix_entry_t entry;
			entry.x = 1;
			i = 0;
			for (entry.j = 0; entry.j < baseMapN2.num; entry.j++)
			{
				// sign is disregarded (bosonic characteristics of pairs)
				if ((baseMapN2.map[entry.j] & baseMapP1.map[n[0]]) != baseMapP1.map[n[0]]) {
					continue;
				}

				// store remaining coordinates after annihilation in tail of 'y'
				FermiDecode(baseMapN2.map[entry.j] - baseMapP1.map[n[0]], y + p2halfTot, N2halfTot - p1halfTot);

				entry.i = Fermi2Base(&baseMapN1, y, N1halfTot);
				if (entry.i == -1) {
					continue;
				}

				assert(i < numEntries);
				matrix[i++] = entry;
			}
			assert(i == numEntries);

			mxSetCell(ret, mxCalcSingleSubscript(ret, 2, n), ToSparseMatrix(baseMapN1.num, baseMapN2.num, matrix, numEntries));

			mxFree(matrix);
		}
	}

	// clean up
	mxFree(y);
	mxFree(p2half);
	mxFree(baseMapP2.map);
	mxFree(baseMapP1.map);
	mxFree(baseMapN2.map);
	mxFree(baseMapN1.map);

	return ret;
}
