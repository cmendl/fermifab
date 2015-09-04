/// \file fermi_map.c
/// \brief Map between bit field representations of fermionic occupations and coordinate lists.
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
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <search.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Encode fermionic coordinates 'x' into a bitfield
///
/// \param x fermionic coordinates
/// \param N length of 'x'
/// \return encoded coordinates
///
bitfield_t FermiEncode(const fermi_coords_t *x, const int N)
{
	int i;
	bitfield_t w = 0;

	for (i = 0; i < N; i++)
	{
		w |= ((bitfield_t)1) << x[i];
	}

	return w;
}

//________________________________________________________________________________________________________________________
///
/// \brief Decode the bitfield 'w' into fermionic coordinates 'x'
///
/// \param w bitfield (input)
/// \param x fermionic coordinates (output)
/// \param N length of 'x'
///
void FermiDecode(const bitfield_t w, fermi_coords_t *x, const int N)
{
	int i = 0, k = 0;

	while (i < N)
	{
		assert(k < 8*sizeof(bitfield_t));
		if (w & (((bitfield_t)1) << k))
		{
			x[i++] = k;
		}
		k++;
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Lexicographically next fermionic bit pattern
///
/// Example: \n
/// f = XXX|01111000_2 -> output XXX|10000111_2,
/// i.e., leading 1 in block of 1s gets shifted to the left,
/// and the remaining 1s are shifted to the end
///
bitfield_t NextFermi(const bitfield_t f)
{
	// compare with 'Bit Twiddling Hacks' by Sean Eron Anderson
	// http://graphics.stanford.edu/~seander/bithacks.html
	bitfield_t t = BitTrailFill(f) + 1;		// -> XXX|10000000_2
	return t | (((LastBit(t) - 1) / LastBit(f)) >> 1);

	// could use 'TZCNT' intrinsic function together with right bit-shift instead of division by LastBit(f)
}


//________________________________________________________________________________________________________________________
///
/// \brief Fermi map bit pattern generation (including configurations);
/// return -1 if last state has been reached
///
/// Encoding: orbital 0 in LSB (least significant bit) \n
///
/// Examples: \n
///
/// orbs = { 8 }, f = XXX|01111000_2 -> output XXX|10000111_2,
/// i.e., leading 1 in block of 1s gets shifted to the left,
/// and the remaining 1s are shifted to the end \n
///
/// orbs = { 6, 5 }, f = 0|01010|111100_2,
/// i.e., particle numbers N = { 4, 2 }
/// -> output 0|01100|001111_2
///
static bitfield_t NextFermiConfig(const int *orbs, const int nc, const bitfield_t f)
{
	assert(orbs[0] > 0);

	bitfield_t t = BitTrailFill(f);		// sets trailing zeros to 1

	// check if 'orbs[0]' orbital group has ended
	bitfield_t mask = (((bitfield_t)1) << orbs[0]) - 1;	// select trailing 'orbs[0]' bits
	if ((t & mask) != mask)		// always false when f == 0
	{
		// compare with 'Bit Twiddling Hacks' by Sean Eron Anderson
		// http://graphics.stanford.edu/~seander/bithacks.html
		t++;	// XXX|10000000_2
		return t | (((LastBit(t) - 1) / LastBit(f)) >> 1);

		// could use 'TZCNT' intrinsic function together with right bit-shift instead of division by LastBit(f)
	}
	else
	{
		if (nc == 1) {
			return (bitfield_t)(-1);
		}

		// store next state of remaining orbital groups in 't'
		t = NextFermiConfig(orbs + 1, nc - 1, f >> orbs[0]);
		if (t == -1) {
			return (bitfield_t)(-1);
		}

		// reset trailing orbital group to |001111_2
		return (mask/LastBit(f)) | (t << orbs[0]);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Create a 'FermiMap': enumerate base indices and map to bit-encoded coordinates
///
/// \param config fermionic configuration
/// \param fm Fermi map (output)
///
void FermiMap(const fermi_config_t *config, fermi_map_t *fm)
{
	int i;
	bitfield_t f = 0;

	fm->num = 1;
	for (i = config->nc - 1; i >= 0; i--)
	{
		assert(0 <= config->N[i] && config->N[i] <= config->orbs[i]);

		// multiply tensor factor dimensions
		fm->num *= Binomial(config->orbs[i], config->N[i]);

		// initialize first coordinates
		f <<= config->orbs[i];
		f |= (((bitfield_t)1) << config->N[i]) - 1;	// init with trailing 'N[i]' bits set to 1
	}
	fm->map = (bitfield_t *)malloc(fm->num*sizeof(bitfield_t));

	// iterate Slater determinants
	for (i = 0; i < fm->num; i++)
	{
		fm->map[i] = f;
		assert(i == 0 || fm->map[i] > fm->map[i-1]);	// should be numerically ordered

		f = NextFermiConfig(config->orbs, config->nc, f);
		assert((f != -1) == (i < fm->num-1));
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Swap two Fermi coordinates
///
static inline void SwapFermiCoords(fermi_coords_t *x, fermi_coords_t *y)
{
	fermi_coords_t t = *x;
	*x = *y;
	*y = t;
}


//________________________________________________________________________________________________________________________
///
/// \brief Sort Fermi coordinates and record permutation sign in { -1, 1 };
/// internal recursive function
///
static void SignedSortInternal(fermi_coords_t *x, const fermi_coords_t left, const fermi_coords_t right, int *sign)
{
	fermi_coords_t i = left, k = right;
	assert(left < right);

	if (x[right] < x[left])
	{
		SwapFermiCoords(&x[left], &x[right]);
		(*sign) = -(*sign);
	}

	while (true)
	{
		while (x[++i] < x[right]) {
			;
		}

		while (x[right] < x[--k]) {
			;
		}

		if (i >= k) {
			break;
		}

		SwapFermiCoords(&x[i], &x[k]);
		(*sign) = -(*sign);
	}

	if (i < right)
	{
		SwapFermiCoords(&x[i], &x[right]);
		(*sign) = -(*sign);
	}

	if (left < i-1)
	{
		SignedSortInternal(x, left, i-1, sign);
	}

	if (i+1 < right)
	{
		SignedSortInternal(x, i+1, right, sign);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Sort coordinates and record permutation sign in { -1, 1 }
///
static inline int SignedSort(fermi_coords_t *x, const int N)
{
	int sign = 1;

	if (N > 1) {
		SignedSortInternal(x, 0, N-1, &sign);
	}

	return sign;
}


//________________________________________________________________________________________________________________________
///
/// \brief Convert Fermi coordinates to base index;
/// return -1 if coordinates are not pairwise different
///
int Fermi2Base(const fermi_map_t *fm, const fermi_coords_t *x, const int N)
{
	// encode 'x'
	bitfield_t w = FermiEncode(x, N);

	if (BitCount(w) < N)
	{
		// elements of 'x' not pairwise different
		return -1;
	}

	bitfield_t *bs = bsearch(&w, fm->map, fm->num, sizeof(bitfield_t), CompareBitfield);
	if (bs == NULL) {
		return -1;
	}

	return (int)(bs - fm->map);
}

//________________________________________________________________________________________________________________________
///
/// \brief Convert Fermi coordinates to base index and record permutation sign;
/// return -1 if coordinates are not pairwise different
///
int Fermi2BaseSign(const fermi_map_t *fm, const fermi_coords_t *x, const int N, int *sign)
{
	int n = Fermi2Base(fm, x, N);

	if (n == -1)
	{
		// no base index match
		(*sign) = 0;
		return -1;
	}

	// obtain permutation sign
	{
		fermi_coords_t *y = (fermi_coords_t *)malloc(N*sizeof(fermi_coords_t));
		memcpy(y, x, N*sizeof(fermi_coords_t));
		(*sign) = SignedSort(y, N);
		free(y);
	}

	return n;
}


//________________________________________________________________________________________________________________________
///
/// \brief Obtain annihilation sign
///
int AnnihilSign(const bitfield_t n, const bitfield_t a)
{
	if ((n & a) != a)
	{
		return 0;
	}
	else
	{
		bitfield_t b = a;	// local copy
		bitfield_t w = n - a;

		// calculate permutation sign
		int count = 0;
		while (b)
		{
			bitfield_t t = LastBit(b);
			// number of particles before particle in current orbital
			count += BitCount(w & (t - 1));
			b -= t;
		}
		return 1 - 2*(count & 1);
	}
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate annihilation sign mask
///
/// Each output bit stores the integer parity of (# following bits) in 'f';
/// example: f = 101100_2 -> output ...11100100_2.
/// Orbital 0 in LSB (least significant bit)
///
bitfield_t AnnihilSignMask(const bitfield_t f)
{
	bitfield_t g = f;	// local copy

	bitfield_t m = 0;
	while (g)
	{
		bitfield_t t = LastBit(g);
		m ^= -t;	// -t: two's complement, e.g. t = 100_2 -> -t = ..111100_2
		g -= t;		// remove last 1-bit from g
	}

	return m;
}
