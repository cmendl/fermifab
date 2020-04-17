/// \file boson_map.c
/// \brief Map between bit field representations of bosonic occupations and coordinate lists.
//
//  Copyright (c) 2008-2015, Christian B. Mendl
//  All rights reserved.
//  http://christian.mendl.net
//
//  This program is free software; you can redistribute it and/or
//  modify it under the terms of the Simplified BSD License
//  http://www.opensource.org/licenses/bsd-license.php
//
//  Reference:
//      Christian B. Mendl
//      The FermiFab toolbox for fermionic many-particle quantum systems
//      Comput. Phys. Commun. 182, 1327-1337 (2011)
//      preprint http://arxiv.org/abs/1103.0872
//________________________________________________________________________________________________________________________
//

#include "boson_map.h"
#include "fermi_map.h"
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <search.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
static int CompareBosonCoords(const void *x, const void *y)
{
	const boson_coords_t u = *(boson_coords_t *)x;
	const boson_coords_t v = *(boson_coords_t *)y;

	if (u < v)
	{
		return -1;
	}
	else if (u == v)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Encode bosonic coordinates 'x' into a bitfield
///
/// Bit pattern 1011001110... means one particle in mode 0, two in mode 1, zero in mode 2, three in mode 3 and zero in mode 4;
/// note that this equals 'FermiEncode(x+{0,1,...N},N)'
///
bitfield_t BosonEncode(const boson_coords_t *x, const int N)
{
	int i;

	// copy input data since we have to sort it
	int *y = (boson_coords_t *)malloc(N*sizeof(boson_coords_t));
	memcpy(y, x, N*sizeof(boson_coords_t));
	qsort(y, N, sizeof(boson_coords_t), CompareBosonCoords);

	bitfield_t w = 0;
	for (i = 0; i < N; i++)
	{
		w |= ((bitfield_t)1) << (y[i] + i);
	}

	free(y);

	return w;
}

//________________________________________________________________________________________________________________________
///
/// \brief Decode the bitfield 'w' into bosonic coordinates 'x'
///
void BosonDecode(const bitfield_t w, boson_coords_t *x, const int N)
{
	int i = 0, k = 0;

	while (i < N)
	{
		assert(k < 8*sizeof(bitfield_t));
		if (w & (((bitfield_t)1) << k))
		{
			x[i] = k - i;
			i++;
		}
		k++;
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Create a 'BosonMap': enumerate base indices and map to bit-encoded coordinates
///
int BosonMap(const int modes, const int N, boson_map_t *bm)
{
	// same as Fermi map - just interpreting bit patterns differently
	int forbs = modes + N - 1;
	fermi_config_t config;
	config.orbs = &forbs;
	config.N    = (int *)&N;
	config.nc   = 1;

	fermi_map_t fm;
	int status = FermiMap(&config, &fm);
	if (status < 0) {
		return status;
	}
	bm->map = fm.map;  // copy pointers
	bm->num = fm.num;

	return 0;
}


//________________________________________________________________________________________________________________________
///
/// \brief Convert bosonic coordinates to base index;
/// return -1 if coordinates cannot be found
///
int Boson2Base(const boson_map_t *bm, const boson_coords_t *x, const int N)
{
	// encode 'x'
	bitfield_t w = BosonEncode(x, N);

	bitfield_t *bs = bsearch(&w, bm->map, bm->num, sizeof(bitfield_t), CompareBitfield);
	if (bs == NULL) {
		return -1;
	}

	return (int)(bs - bm->map);
}
