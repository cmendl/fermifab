/// \file gen_comprise.c
/// \brief Create index mapping from input to target configuration; Matlab 'mex' function.
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
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief X = gen_comprise(orbs,N,merge_next)
///
/// Comprise input configuration scope, i.e., create index mapping from input to target configuration
/// defined by merging subspace neighbors in the input configuration. For example, the input
/// orbs = {3,4,6}, N = {1,2,4} with merge_next = {true,false} results in orbs = {7,6}, N = {3,4};
/// Matlab 'mex' function
///
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	int i;

	if (nlhs > 1)  { mexErrMsgTxt("One output argument expected."); }
	if (nrhs != 3) { mexErrMsgTxt("Three input arguments required."); }

	// number of orbitals
	fermi_config_t config;
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

	// particle configuration
	if ((mxGetM(prhs[1])!=config.nc || mxGetN(prhs[1])!=1) && (mxGetM(prhs[1])!=1 || mxGetN(prhs[1])!=config.nc)) {
		mexErrMsgTxt("Number of orbitals and particle fields must have the same length.");
	}
	config.N = (int *)mxMalloc(config.nc*sizeof(int));
	for (i = 0; i < config.nc; i++)
	{
		config.N[i] = (int)(*(mxGetPr(prhs[1])+i));
		if (config.N[i] < 0 || config.N[i] > config.orbs[i]) {
			mexErrMsgTxt("Particle numbers must be between zero and the corresponding number of orbitals.");
		}
	}

	// read merge table
	if ((mxGetM(prhs[2])!=config.nc-1 || mxGetN(prhs[2])!=1) && (mxGetM(prhs[2])!=1 || mxGetN(prhs[2])!=config.nc-1)) {
		mexErrMsgTxt("Merge table must have length of 'orbs' field - 1.");
	}
	bool *merge_next = (bool *)mxMalloc((config.nc-1)*sizeof(bool));
	if (mxIsLogical(prhs[2]))
	{
		for (i = 0; i < config.nc-1; i++)
		{
			merge_next[i] = *(mxGetLogicals(prhs[2])+i);
		}
	}
	else
	{
		for (i = 0; i < config.nc-1; i++)
		{
			merge_next[i] = (*(mxGetPr(prhs[2])+i) != 0);
		}
	}

	// generate 'merged' orbs and particle numbers
	fermi_config_t config_merged;
	config_merged.orbs = (int *)mxMalloc(config.nc*sizeof(int));
	config_merged.N    = (int *)mxMalloc(config.nc*sizeof(int));
	config_merged.orbs[0] = config.orbs[0];
	config_merged.N   [0] = config.N   [0];
	config_merged.nc = 1;
	for (i = 0; i < config.nc-1; i++)
	{
		if (merge_next[i])
		{
			config_merged.orbs[config_merged.nc-1] += config.orbs[i+1];
			config_merged.N   [config_merged.nc-1] += config.N   [i+1];
		}
		else
		{
			config_merged.nc++;
			config_merged.orbs[config_merged.nc-1] = config.orbs[i+1];
			config_merged.N   [config_merged.nc-1] = config.N   [i+1];
		}
	}

	fermi_map_t baseMap, baseMap_merged;
	FermiMap(&config, &baseMap);
	FermiMap(&config_merged, &baseMap_merged);

	// output data setup
	// allocate sparse matrix
	plhs[0] = mxCreateSparse(baseMap_merged.num, baseMap.num, baseMap.num, mxREAL);
	// fill 'jc'
	mwIndex *jc = mxGetJc(plhs[0]);
	for (i = 0; i <= baseMap.num; i++)
	{
		jc[i] = i;
	}
	// set indices
	double  *sr = mxGetPr(plhs[0]);
	mwIndex *ir = mxGetIr(plhs[0]);
	for (i = 0; i < baseMap.num; i++)
	{
		bitfield_t *bs = bsearch(&baseMap.map[i], baseMap_merged.map, baseMap_merged.num, sizeof(bitfield_t), CompareBitfield);
		assert(bs != NULL);

		ir[i] = bs - baseMap_merged.map;
		sr[i] = 1;
	}

	// clean up
	mxFree(baseMap_merged.map);
	mxFree(config_merged.N);
	mxFree(config_merged.orbs);
	mxFree(merge_next);
	mxFree(baseMap.map);
	mxFree(config.N);
	mxFree(config.orbs);
}
