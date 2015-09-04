/// \file boson2coords.c
/// \brief Convert bosonic basis index to lexicographically ordered N-particle 'coordinates'; Matlab 'mex' function.
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
#include "boson_map.h"
#include <mex.h>


//________________________________________________________________________________________________________________________
///
/// \brief X = boson2coords(modes,N)
///
/// Convert bosonic basis index to lexicographically ordered N-particle 'coordinates';
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nlhs > 1)  { mexErrMsgTxt("One output argument expected."); }
	if (nrhs != 2) { mexErrMsgTxt("Two input arguments required."); }

	int modes = (int)mxGetScalar(prhs[0]);
	int N     = (int)mxGetScalar(prhs[1]);

	if (modes < 1) { mexErrMsgTxt("Number of 'modes' must at least be 1."); }
	if (N < 1)     { mexErrMsgTxt("'N' must be greater than zero."); }

	boson_map_t baseMap;
	BosonMap(modes, N, &baseMap);
	boson_coords_t *x = (boson_coords_t *)mxMalloc(N*sizeof(boson_coords_t));

	// copy output data
	plhs[0] = mxCreateDoubleMatrix(N, baseMap.num, mxREAL);
	double *d = mxGetPr(plhs[0]);
	mwIndex subs[2];
	for (subs[1] = 0; subs[1] < (mwIndex)baseMap.num; subs[1]++)
	{
		// get coordinates
		BosonDecode(baseMap.map[subs[1]], x, N);

		for (subs[0] = 0; subs[0] < (mwIndex)N; subs[0]++) {
			d[mxCalcSingleSubscript(plhs[0], 2, subs)] = (double)(x[subs[0]]+1);	// convert between C++ zero-based indices and Matlab 1-based indices
		}
	}

	// clean up
	mxFree(baseMap.map);
	mxFree(x);
}
