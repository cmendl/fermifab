/// \file fermi2coords.c
/// \brief Convert Slater basis index to lexicographically ordered N-particle orbital 'coordinates'; Matlab 'mex' function.
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
#include <mex.h>


//________________________________________________________________________________________________________________________
///
/// \brief X = fermi2coords(orbs,N)
///
/// Convert Slater basis index to lexicographically ordered N-particle orbital 'coordinates';
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;

	if (nlhs > 1)  { mexErrMsgTxt("One output argument expected."); }
	if (nrhs != 2) { mexErrMsgTxt("Two input arguments required."); }

	fermi_config_t config;

	// number of orbitals
	if (mxGetM(prhs[0]) == 1)
	{
		config.nc = (int)mxGetN(prhs[0]);
	}
	else
	{
		config.nc = (int)mxGetM(prhs[0]);
		if (mxGetN(prhs[0]) != 1) {
			mexErrMsgTxt("Invalid 'orbs' input argument.");
		}
	}
	if (config.nc == 0) { mexErrMsgTxt("No 'orbs' input argument."); }
	config.orbs = (int *)mxMalloc(config.nc*sizeof(int));
	for (i = 0; i < config.nc; i++)
	{
		config.orbs[i] = (int)(*(mxGetPr(prhs[0])+i));

		if (config.orbs[i] < 1) {
			mexErrMsgTxt("Number of orbitals must be positive.");
		}
	}

	// particle numbers
	if ((mxGetM(prhs[1])!=config.nc || mxGetN(prhs[1])!=1) && (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=config.nc)) {
		mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
	}
	config.N = (int *)mxMalloc(config.nc*sizeof(int));
	for (i = 0; i < config.nc; i++)
	{
		config.N[i] = (int)(*(mxGetPr(prhs[1])+i));

		if (config.N[i] < 0 || config.N[i] > config.orbs[i]) {
			mexErrMsgTxt("Particle numbers must be between 0 and the corresponding number of orbitals.");
		}
	}
	int Ntot = IntegerSum(config.N, config.nc);

	fermi_map_t baseMap;
	FermiMap(&config, &baseMap);
	fermi_coords_t *x = (fermi_coords_t *)mxMalloc(Ntot*sizeof(fermi_coords_t));

	// copy output data
	plhs[0] = mxCreateDoubleMatrix(Ntot, baseMap.num, mxREAL);
	double *d = mxGetPr(plhs[0]);
	mwIndex subs[2];
	for (subs[1] = 0; subs[1] < (mwIndex)baseMap.num; subs[1]++)
	{
		// get coordinates
		FermiDecode(baseMap.map[subs[1]], x, Ntot);

		for (subs[0] = 0; subs[0] < (mwIndex)Ntot; subs[0]++)
		{
			d[mxCalcSingleSubscript(plhs[0], 2, subs)] = (double)(x[subs[0]] + 1);	// convert between C++ zero-based indices and Matlab 1-based indices
		}
	}

	// clean up
	mxFree(x);
	mxFree(baseMap.map);
	mxFree(config.N);
	mxFree(config.orbs);
}
