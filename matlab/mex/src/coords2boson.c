/// \file coords2boson.c
/// \brief Convert lexicographically ordered N-particle 'coordinates' to the bosonic basis index; Matlab 'mex' function.
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

#include "boson_map.h"
#include "util.h"
#include <mex.h>


//________________________________________________________________________________________________________________________
///
/// \brief m = coords2boson(x,modes)
///
/// Convert N-particle 'coordinates' to the corresponding bosonic basis index;
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;

	if (nrhs != 2) { mexErrMsgTxt("2 input arguments required."); }
	if (nlhs > 2)  { mexErrMsgTxt("Wrong number of output arguments."); }

	// number of particles
	int N;
	if (mxGetM(prhs[0]) == 1)
	{
		N = (int)mxGetN(prhs[0]);
	}
	else
	{
		N = (int)mxGetM(prhs[0]);
		if (mxGetN(prhs[0]) != 1) {
			mexErrMsgTxt("'x' must be a vector of integers.");
		}
	}
	if (N == 0) { mexErrMsgTxt("Invalid 'x' input argument."); }

	int modes = (int)mxGetScalar(prhs[1]);
	if (modes < 1) { mexErrMsgTxt("Number of 'modes' must at least be 1."); }

	boson_coords_t *x = (boson_coords_t *)mxMalloc(N*sizeof(boson_coords_t));
	double *d = mxGetPr(prhs[0]);
	for (i = 0; i < N; i++)
	{
		x[i] = (int)d[i]-1;	// convert between C++ zero-based indices and Matlab 1-based indices

		if (x[i] < 0 || x[i] >= modes) {
			mexErrMsgTxt("Invalid input argument 'x': out of range.");
		}
	}

	// calculate base index
	boson_map_t baseMap;
	BosonMap(modes, N, &baseMap);

	int m = Boson2Base(&baseMap, x, N);
	// indices start at 1 in Matlab
	if (m != -1) {
		m++;
	}

	// copy output data
	plhs[0] = mxCreateDoubleScalar(m);

	// clean up
	mxFree(baseMap.map);
	mxFree(x);
}
