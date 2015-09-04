/// \file gen_spintrace_coul.c
/// \brief Trace out spin from 2-particle Coulomb operator for '2*spatialorbs' orbitals including spin,
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

#include "util.h"
#include "mextools.h"
#include "fermi_map.h"
#include "boson_map.h"
#include <assert.h>


static mxArray *GenerateSpintraceCoulomb(const int spatialorbs);


//________________________________________________________________________________________________________________________
///
/// \brief K = gen_spintrace_coul(spatialorbs)
///
/// Trace out spin from 2-particle Coulomb operator; Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs > 1)  { mexErrMsgTxt("1 output argument expected."); }
	if (nrhs != 1) { mexErrMsgTxt("1 input arguments required."); }

	int spatialorbs = (int)mxGetScalar(prhs[0]);
	if (spatialorbs < 1) { mexErrMsgTxt("Number of spatial orbitals 'spatialorbs' must at least be 1."); }

	// verbose output
	// mexPrintf("Calculating spin trace Coulomb mapping for %i orbitals (including spin)...", 2*spatialorbs);
	plhs[0] = GenerateSpintraceCoulomb(spatialorbs);
	// mexPrintf(" Done.\n");
}


//________________________________________________________________________________________________________________________
///
/// \brief Swap two bosonic coordinates
///
static inline void SwapBosonCoords(boson_coords_t *x, boson_coords_t *y)
{
	boson_coords_t t = *x;
	*x = *y;
	*y = t;
}


//________________________________________________________________________________________________________________________
///
/// \brief Apply symmetry to Coulomb integrals (ab|cd) with a <= b, c <= d, (a,b) <= (c,d) according to Bose map ordering
///
static void ApplySymmetry(boson_coords_t *ab, boson_coords_t *cd)
{
	// assume that all spatial orbitals are real,
	// so e.g. (ab|cd) = (ba|cd)
	if (ab[0] > ab[1]) { SwapBosonCoords(&ab[0], &ab[1]); }
	if (cd[0] > cd[1]) { SwapBosonCoords(&cd[0], &cd[1]); }

	// Coulomb operator is symmetric;
	// use reverse-lexicographical Bose map ordering (last entry more significant than first)
	if (ab[1] > cd[1])
	{
		SwapBosonCoords(&ab[0], &cd[0]);
		SwapBosonCoords(&ab[1], &cd[1]);
	}
	else if (ab[1] == cd[1] && ab[0] > cd[0])
	{
		SwapBosonCoords(&ab[0], &cd[0]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Trace out spin from 2-particle Coulomb operator for '2*spatialorbs' orbitals including spin,
/// i.e. the orbitals are (1 up, 1 down, 2 up, 2 down,...) with spatial orbitals 1, 2,...
///
static mxArray *GenerateSpintraceCoulomb(const int spatialorbs)
{
	int i;
	mwIndex m[2];

	// create Bose map
	boson_map_t coulombMap;
	BosonMap(spatialorbs, 2, &coulombMap);
	// create Fermi map
	fermi_map_t baseMap2;
	{
		int orbs = 2*spatialorbs;
		int fN = 2;

		fermi_config_t config;
		config.orbs = &orbs;
		config.N    = &fN;
		config.nc   = 1;
		FermiMap(&config, &baseMap2);
	}

	// allocate cell matrix array
	mxArray *ret = mxCreateCellMatrix(coulombMap.num, coulombMap.num);
	if (!ret) { mexErrMsgTxt("Creating cell matrix failed."); }

	// fill cells with our own 'sparse' matrix structure first
	sparse_matrix_entry_t **cells = (sparse_matrix_entry_t **)mxMalloc(coulombMap.num*coulombMap.num*sizeof(sparse_matrix_entry_t *));
	int *cellNum                  =                    (int *)mxCalloc(coulombMap.num*coulombMap.num, sizeof(int));

	sparse_matrix_entry_t entry;
	for (entry.j = 0; entry.j < baseMap2.num; entry.j++)
	{
		// fermionic coordinates
		fermi_coords_t p[2];
		FermiDecode(baseMap2.map[entry.j], p, 2);

		for (entry.i = 0; entry.i < baseMap2.num; entry.i++)
		{
			// fermionic coordinates
			fermi_coords_t q[2];
			FermiDecode(baseMap2.map[entry.i], q, 2);

			// test if both spins match
			if (((p[0] - q[0]) % 2 == 0) && ((p[1] - q[1]) % 2 == 0))
			{
				boson_coords_t ab[2], cd[2];	// spatial orbitals in Coulomb integrals

				// (ab|cd)
				ab[0] = (boson_coords_t)(p[0]/2);	// round down
				ab[1] = (boson_coords_t)(q[0]/2);
				cd[0] = (boson_coords_t)(p[1]/2);
				cd[1] = (boson_coords_t)(q[1]/2);
				ApplySymmetry(ab, cd);

				// positive sign
				entry.x = 1;

				// map to bosonic indices of return value
				m[0] = Boson2Base(&coulombMap, ab, 2);
				m[1] = Boson2Base(&coulombMap, cd, 2);
				i = (int)mxCalcSingleSubscript(ret, 2, m);
				if (cellNum[i] == 0)
				{
					cells[i] = (sparse_matrix_entry_t *)mxMalloc(4*sizeof(sparse_matrix_entry_t));
				}
				else if (cellNum[i] == 4)
				{
					cells[i] = (sparse_matrix_entry_t *)mxRealloc(cells[i], 16*sizeof(sparse_matrix_entry_t));
				}
				cells[i][cellNum[i]++] = entry;
			}

			if (((p[0] - q[1]) % 2 == 0) && ((p[1] - q[0]) % 2 == 0))
			{
				boson_coords_t ab[2], cd[2];	// spatial orbitals in Coulomb integrals

				// (ab|cd)
				ab[0] = (boson_coords_t)(p[0]/2);	// round down
				ab[1] = (boson_coords_t)(q[1]/2);
				cd[0] = (boson_coords_t)(p[1]/2);
				cd[1] = (boson_coords_t)(q[0]/2);
				ApplySymmetry(ab, cd);

				// minus sign in exchange integral
				entry.x = -1;

				// map to bosonic indices of return value
				m[0] = Boson2Base(&coulombMap, ab, 2);
				m[1] = Boson2Base(&coulombMap, cd, 2);
				i = (int)mxCalcSingleSubscript(ret, 2, m);
				if (cellNum[i] == 0)
				{
					cells[i] = (sparse_matrix_entry_t *)mxMalloc(4*sizeof(sparse_matrix_entry_t));
				}
				else if (cellNum[i] == 4)
				{
					cells[i] = (sparse_matrix_entry_t *)mxRealloc(cells[i], 16*sizeof(sparse_matrix_entry_t));
				}
				cells[i][cellNum[i]++] = entry;
			}
		}
	}

	// create Matlab sparse structure
	for (m[0] = 0; m[0] < (mwIndex)coulombMap.num; m[0]++)
	{
		for (m[1] = 0; m[1] < (mwIndex)coulombMap.num; m[1]++)
		{
			i = (int)mxCalcSingleSubscript(ret, 2, m);
			assert(cellNum[i] < 17);
			mxSetCell(ret, i, ToSparseMatrix(baseMap2.num, baseMap2.num, cells[i], cellNum[i]));
		}
	}

	// tidy up
	for (i = 0; i < coulombMap.num*coulombMap.num; i++)
	{
		if (cellNum[i] > 0) {
			mxFree(cells[i]);
		}
	}
	mxFree(cells);
	mxFree(cellNum);
	mxFree(baseMap2.map);
	mxFree(coulombMap.map);

	return ret;
}
