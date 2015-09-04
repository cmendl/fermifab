/// \file coords2fermi.c
/// \brief Convert lexicographically ordered N-particle orbital 'coordinates' to a Slater basis index; Matlab 'mex' function.
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
#include <memory.h>
#include <mex.h>


//________________________________________________________________________________________________________________________
///
/// \brief n = coords2fermi(x,orbs) \n
/// [n,sign] = coords2fermi(x,orbs)
///
/// Convert N-particle orbital 'coordinates' to the corresponding Slater basis index;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i, j, m;

	if (nrhs != 2) { mexErrMsgTxt("2 input arguments required."); }
	if (nlhs > 2)  { mexErrMsgTxt("Wrong number of output arguments."); }

	// total number of particles
	int Ntot;
	if (mxGetM(prhs[0]) == 1)
	{
		Ntot = (int)mxGetN(prhs[0]);
	}
	else
	{
		Ntot = (int)mxGetM(prhs[0]);
		if (mxGetN(prhs[0]) != 1) {
			mexErrMsgTxt("'x' must be a vector of integers.");
		}
	}
	if (Ntot == 0) { mexErrMsgTxt("Invalid 'x' input argument."); }

	fermi_config_t config;

	// number of orbitals
	if (mxGetM(prhs[1]) == 1)
	{
		config.nc = (int)mxGetN(prhs[1]);
	}
	else
	{
		config.nc = (int)mxGetM(prhs[1]);
		if (mxGetN(prhs[1]) != 1) {
			mexErrMsgTxt("Invalid 'orbs' input argument.");
		}
	}
	if (config.nc == 0) { mexErrMsgTxt("No 'orbs' input argument."); }
	config.orbs = (int *)mxMalloc(config.nc*sizeof(int));
	for (i = 0; i < config.nc; i++)
	{
		config.orbs[i] = (int)(*(mxGetPr(prhs[1])+i));

		if (config.orbs[i] < 1) {
			mexErrMsgTxt("Number of orbitals must be positive.");
		}
	}

	// number of particles in the individual orbital groups
	config.N = (int *)mxMalloc(config.nc*sizeof(int));
	memset(config.N, 0, config.nc*sizeof(int));

	fermi_coords_t *x = (fermi_coords_t *)mxMalloc(Ntot*sizeof(fermi_coords_t));
	for (i = 0; i < Ntot; i++)
	{
		x[i] = (int)(*(mxGetPr(prhs[0])+i))-1;	// convert between C++ zero-based indices and Matlab 1-based indices

		if (x[i] < 0) {
			mexErrMsgTxt("Input argument 'x' must only contain positive entries.");
		}

		// increase the number of particles in corresponding orbital group
		m = 0;
		for (j = 0; j < config.nc; j++)
		{
			m += config.orbs[j];
			if (x[i] < m)
			{
				config.N[j]++;
				break;
			}
		}
		if (j == config.nc) {
			mexErrMsgTxt("Input argument 'x' out of range.");
		}
	}

	for (i = 0; i < config.nc; i++)
	{
		if (config.N[i] > config.orbs[i]) {
			mexErrMsgTxt("Number of orbitals must at least be equal to the number of particles.");
		}
	}

	// build mapping table
	fermi_map_t baseMap;
	FermiMap(&config, &baseMap);

	if (nlhs == 2)
	{
		// base index and sign
		int sign;
		m = Fermi2BaseSign(&baseMap, x, Ntot, &sign);

		plhs[1] = mxCreateDoubleScalar(sign);
	}
	else
	{
		// base index
		m = Fermi2Base(&baseMap, x, Ntot);
	}

	// copy base index
	if (m != -1) {
		m++;	// indices start at 1 in Matlab
	}
	plhs[0] = mxCreateDoubleScalar(m);

	// clean up
	mxFree(baseMap.map);
	mxFree(x);
	mxFree(config.N);
	mxFree(config.orbs);
}
