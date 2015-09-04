/// \file map_pairs.c
/// \brief Map pair structure Slater basis indices to full structure indices; Matlab 'mex' function.
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


//________________________________________________________________________________________________________________________
///
/// \brief P = map_pairs(orbshalf,Nhalf)
///
/// Map pair structure Slater basis indices to full structure indices;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;

	if (nrhs != 2) { mexErrMsgTxt("Two input arguments required."); }
	if (nlhs > 1)  { mexErrMsgTxt("Wrong number of output arguments."); }

	int orbshalf = (int)mxGetScalar(prhs[0]);
	int Nhalf    = (int)mxGetScalar(prhs[1]);

	if (Nhalf < 1)        { mexErrMsgTxt("'Nhalf' must at least be 1."); }
	if (Nhalf > orbshalf) { mexErrMsgTxt("'orbshalf' must at least be 'Nhalf'."); }

	int orbs = 2*orbshalf;
	int N    = 2*Nhalf;

	// create 'half' Fermi map
	fermi_map_t baseMapHalf;
	{
		fermi_config_t config;
		config.orbs = &orbshalf;
		config.N    = &Nhalf;
		config.nc   = 1;
		FermiMap(&config, &baseMapHalf);
	}

	// create 'full' Fermi map
	fermi_map_t baseMap;
	{
		fermi_config_t config;
		config.orbs = &orbs;
		config.N    = &N;
		config.nc   = 1;
		FermiMap(&config, &baseMap);
	}

	// note that P is an isometry
	sparse_matrix_entry_t *P = (sparse_matrix_entry_t *)mxMalloc(baseMapHalf.num*sizeof(sparse_matrix_entry_t));
	sparse_matrix_entry_t entry;
	entry.x = 1;

	fermi_coords_t *x     = (fermi_coords_t *)mxMalloc(N    *sizeof(fermi_coords_t));
	fermi_coords_t *xhalf = (fermi_coords_t *)mxMalloc(Nhalf*sizeof(fermi_coords_t));

	// map base indices
	for (entry.j = 0; entry.j < baseMapHalf.num; entry.j++)
	{
		FermiDecode(baseMapHalf.map[entry.j], xhalf, Nhalf);
		for (i = 0; i < Nhalf; i++)
		{
			x[2*i  ] = 2*xhalf[i];
			x[2*i+1] = 2*xhalf[i]+1;
		}
		entry.i = Fermi2Base(&baseMap, x, Nhalf);

		P[entry.j] = entry;
	}

	// copy output data
	plhs[0] = ToSparseMatrix(baseMap.num, baseMapHalf.num, P, baseMapHalf.num);

	mxFree(baseMap.map);
	mxFree(baseMapHalf.map);
	mxFree(xhalf);
	mxFree(x);
}
