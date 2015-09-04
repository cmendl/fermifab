/// \file gen_rdm.c
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices; Matlab 'mex' function.
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

#include "util.h"
#include "fermi_map.h"
#include "mextools.h"
#include <assert.h>


static mxArray *GenerateRDM(const int *orbs, const int *p1, const int *N1, const int *N2, const int nc);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_rdm(orbs,p,N) \n
///        K = gen_rdm(orbs,p1,N1,N2)
///
/// Create the matrix kernels K required for calculating p-body reduced density matrices;
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
			mexErrMsgTxt("Invalid 'orbs' input argument.");
		}
	}
	if (nc == 0) { mexErrMsgTxt("No 'orbs' input argument."); }
	int *orbs = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		orbs[i] = (int)(*(mxGetPr(prhs[0])+i));
		if (orbs[i] < 1) {
			mexErrMsgTxt("Number of orbitals must be positive.");
		}
	}

	// reduced particle configuration
	if ((mxGetM(prhs[1])!=nc || mxGetN(prhs[1])!=1) && (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=nc)) {
		mexErrMsgTxt("Number of orbitals and reduced particle fields must have the same length.");
	}
	int ptot = 0;	// total reduced particle number
	int *p1 = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p1[i] = (int)(*(mxGetPr(prhs[1])+i));
		if (p1[i] < 0) { mexErrMsgTxt("No negative particle numbers allowed."); }
		ptot += p1[i];
	}
	if (ptot < 1) { mexErrMsgTxt("Reduced particle number must at least be 1."); }

	// particle configurations of full fermionic spaces
	if ((mxGetM(prhs[2])!=nc || mxGetN(prhs[2])!=1) && (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=nc)) {
		mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
	}
	int Ntot = 0;	// total particle number
	int *N1 = (int *)mxMalloc(nc*sizeof(int));
	int *N2;
	for (i = 0; i < nc; i++)
	{
		N1[i] = (int)(*(mxGetPr(prhs[2])+i));
		if (N1[i] < 0 || N1[i] > orbs[i]) {
			mexErrMsgTxt("Particle numbers must be between zero and the corresponding number of orbitals.");
		}
		Ntot += N1[i];
	}
	if (nrhs == 4)
	{
		if ((mxGetM(prhs[3])!=nc || mxGetN(prhs[3])!=1) && (mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=nc)) {
			mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
		}

		int Ndiff = 0;
		N2 = (int *)mxMalloc(nc*sizeof(int));
		for (i = 0; i < nc; i++)
		{
			N2[i] = (int)(*(mxGetPr(prhs[3])+i));
			if (N2[i] < 0 || N2[i] > orbs[i]) {
				mexErrMsgTxt("Particle numbers must be between zero and the corresponding number of orbitals.");
			}
			Ndiff += N1[i] - N2[i];
		}
		// check whether total particle numbers agree
		if (Ndiff != 0) {
			mexErrMsgTxt("Total particle numbers in 'N1' and 'N2' must be the same.");
		}
	}
	else
	{
		N2 = N1;
	}

	// consistency checks
	if (ptot > Ntot) { mexErrMsgTxt("Reduced particle number must be less or equal to the number of particles in the full system."); }
	for (i = 0; i < nc; i++)
	{
		if (p1[i] > N2[i]) {
			mexErrMsgTxt("'p1' must be less or equal to the number of particles in the full system.");
		}
		if (N1[i] - N2[i] + p1[i] < 0)	// p2[i]
		{
			mexErrMsgTxt("No negative particle numbers allowed for p2 := N1 - N2 + p1.");
		}
	}

	// verbose output
	// mexPrintf("Creating reduced %i-body matrices from %i-fermion states (number of orbitals: %i)...", ptot, Ntot, IntegerSum(orbs, nc));
	plhs[0] = GenerateRDM(orbs, p1, N1, N2, nc);
	// mexPrintf(" Done.\n");

	// clean up
	mxFree(N1);
	if (N1 != N2) {
		mxFree(N2);
	}
	mxFree(p1);
	mxFree(orbs);
}


//________________________________________________________________________________________________________________________
///
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices
/// on a Hilbert space with dimension 'orbs'
///
/// Given a normalized N-body wavefunction psi (coefficients of canonical Slater basis in lexicographical ordering),
/// <psi | K{i,j} psi> is the coefficient i,j of the p-body reduced density matrix of psi.
/// In the creation/annihilation formalism, K{i,j} is the operator a_j^dagger a_i acting on the N-body space wedge^N H.
///
static mxArray *GenerateRDM(const int *orbs, const int *p1, const int *N1, const int *N2, const int nc)
{
	int i;

	const int N1tot = IntegerSum(N1, nc);
	const int N2tot = IntegerSum(N2, nc);
	const int p1tot = IntegerSum(p1, nc);
	const int p2tot =  N1tot - N2tot + p1tot;

	// basic configuration setup
	fermi_config_t config;
	config.orbs = (int *)orbs;
	config.nc = nc;

	// fermionic particle spaces with configurations N1 and N2
	fermi_map_t baseMapN1, baseMapN2;
	config.N = (int *)N1; FermiMap(&config, &baseMapN1);
	config.N = (int *)N2; FermiMap(&config, &baseMapN2);

	// configurations p1 and p2
	fermi_map_t baseMapP1, baseMapP2;
	config.N = (int *)p1; FermiMap(&config, &baseMapP1);
	// N2 - p1 == N1 - p2
	int *p2 = (int *)mxMalloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p2[i] = N1[i] - N2[i] + p1[i];
		assert(p2[i] >= 0);
	}
	config.N = p2; FermiMap(&config, &baseMapP2);

	// allocate cell matrix array
	mxArray *ret = mxCreateCellMatrix(baseMapP1.num, baseMapP2.num);
	if (!ret) { mexErrMsgTxt("Creating cell matrix failed (probably out of memory)."); }

	fermi_coords_t *y = (fermi_coords_t *)mxMalloc(N1tot*sizeof(fermi_coords_t));

	mwIndex n[2];
	for (n[0] = 0; n[0] < (mwIndex)baseMapP1.num; n[0]++)
	{
		for (n[1] = 0; n[1] < (mwIndex)baseMapP2.num; n[1]++)
		{
			// fill first 'p2tot' coordinates of 'y' only
			FermiDecode(baseMapP2.map[n[1]], y, p2tot);

			// determine number of matrix entries for current p1 and p2 configuration
			int numEntries = 1;
			bitfield_t conf = baseMapP1.map[n[0]] | baseMapP2.map[n[1]];
			for (i = 0; i < nc; i++)
			{
				bitfield_t mask = ((bitfield_t)1 << orbs[i]) - 1;
				numEntries *= Binomial(orbs[i] - BitCount(mask & conf), N2[i] - p1[i]);
				conf <<= orbs[i];
			}

			sparse_matrix_entry_t *matrix = (sparse_matrix_entry_t *)mxMalloc(numEntries*sizeof(sparse_matrix_entry_t));

			i = 0;
			sparse_matrix_entry_t entry;
			for (entry.j = 0; entry.j < baseMapN2.num; entry.j++)
			{
				int sign[2];
				sign[0] = AnnihilSign(baseMapN2.map[entry.j], baseMapP1.map[n[0]]);
				if (!sign[0]) {
					continue;
				}

				// store remaining coordinates after annihilation in tail of 'y'
				FermiDecode(baseMapN2.map[entry.j] - baseMapP1.map[n[0]], y + p2tot, N2tot - p1tot);

				entry.i = Fermi2BaseSign(&baseMapN1, y, N1tot, &sign[1]);
				if (entry.i == -1) {
					continue;
				}

				entry.x = sign[0]*sign[1];

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
	mxFree(p2);
	mxFree(baseMapP2.map);
	mxFree(baseMapP1.map);
	mxFree(baseMapN2.map);
	mxFree(baseMapN1.map);

	return ret;
}
