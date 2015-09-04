/// \file gen_config.c
/// \brief Enumerate configurations; Matlab 'mex' function.
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
#include <memory.h>
#include <assert.h>


static bool NextConfig(const int nc, int *conf);
static bool NextConfigMax(const int nc, const int *m, int *conf);


//________________________________________________________________________________________________________________________
///
/// \brief C = gen_config(nc,N) \n
///        C = gen_config(nc,N,maxC)
///
/// Enumerate all configurations of 'N' particles in 'nc' slots,
/// such that not more than 'maxC[i]' particles are in slot 'i';
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;
	int nc, N, left;
	double *d;
	int *conf, *confnum, *maxC;
	mwIndex subs[2];

	if (nlhs > 1) { mexErrMsgTxt("One output argument expected."); }
	if (nrhs < 2 || nrhs > 3) { mexErrMsgTxt("Two or three input arguments required."); }

	// number of slots
	nc = (int)mxGetScalar(prhs[0]);
	if (nc < 1) { mexErrMsgTxt("Number of slots 'nc' must at least be 1."); }

	// particle number
	N = (int)mxGetScalar(prhs[1]);
	if (N < 0) { mexErrMsgTxt("'N' must be non-negative."); }

	conf = (int *)mxCalloc(nc, sizeof(int));

	if (nrhs == 2)
	{
		plhs[0] = mxCreateDoubleMatrix(nc, Binomial(nc + N - 1, N), mxREAL);
		d = mxGetPr(plhs[0]);

		conf[0] = N;
		subs[1] = 0;
		do
		{
			for (subs[0] = 0; subs[0] < (mwIndex)nc; subs[0]++)
			{
				d[mxCalcSingleSubscript(plhs[0], 2, subs)] = conf[subs[0]];
			}

			subs[1]++;
		}
		while (NextConfig(nc, conf));

		assert(subs[1] == mxGetN(plhs[0]));
	}
	else
	{
		// get maxmum number of particles in each slot
		if ((mxGetM(prhs[2])!=nc || mxGetN(prhs[2])!=1) && (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=nc)) {
			mexErrMsgTxt("Input argument 'maxC' must have length 'nc'.");
		}
		maxC = (int *)mxMalloc(nc*sizeof(int));
		for (i = 0; i < nc; i++)
		{
			maxC[i] = (int)(*(mxGetPr(prhs[2])+i));
			if (maxC[i] < 0) {
				mexErrMsgTxt("Maximum number of particles in each slot must be non-negative.");
			}
		}

		// first valid configuration
		left = N;
		for (i = 0; i < nc; i++)
		{
			conf[i] = mini(left, maxC[i]);
			left -= conf[i];
		}
		if (left > 0)	// no single valid configuration
		{
			plhs[0] = mxCreateDoubleMatrix(nc, 0, mxREAL);
			mxFree(maxC);
			mxFree(conf);
			return;
		}

		// determine total number of valid configurations
		confnum = (int *)mxMalloc(nc*sizeof(int));
		memcpy(confnum, conf, nc*sizeof(int));
		i = 1;
		while (NextConfigMax(nc, maxC, confnum))
		{
			i++;
		}
		mxFree(confnum);

		// copy configurations to output array
		plhs[0] = mxCreateDoubleMatrix(nc, i, mxREAL);
		d = mxGetPr(plhs[0]);
		subs[1] = 0;
		do
		{
			for (subs[0] = 0; subs[0] < (mwIndex)nc; subs[0]++)
			{
				d[mxCalcSingleSubscript(plhs[0], 2, subs)] = conf[subs[0]];
			}

			subs[1]++;
		}
		while (NextConfigMax(nc, maxC, conf));

		mxFree(maxC);
	}

	mxFree(conf);
}


//________________________________________________________________________________________________________________________
///
/// \brief Generate next configuration
///
static bool NextConfig(const int nc, int *conf)
{
	int i;
	assert(nc > 0);

	if (nc == 1) {
		return false;
	}

	if (conf[0] > 0)
	{
		conf[0]--;
		conf[1]++;

		return true;
	}
	else
	{
		i = 1;
		while (conf[i] == 0)
		{
			i++;
			assert(i < nc);
		}
		conf[0] = conf[i] - 1;
		conf[i] = 1;

		return NextConfig(nc - i, conf + i);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Generate next configuration given the upper limits in 'm'
/// by filtering valid configurations
///
static bool NextConfigMax(const int nc, const int *m, int *conf)
{
	int i;
	do
	{
		if (!NextConfig(nc, conf)) {
			return false;
		}

		// check if configuration is valid, i.e.,
		// whether conf[i] <= m[i] for all i
		i = 0;
		do
		{
			if (conf[i] > m[i]) {
				break;
			}
			i++;
		}
		while (i < nc);
	}
	while (i < nc);

	return true;
}
