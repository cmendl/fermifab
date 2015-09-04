/// \file spintrace_coul.c
/// \brief Trace out spin and determine to-be calculated Coulomb integrals;
/// use fermionic coordinates as input, and encode spatial orbitals as bosonic bitfields.
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

#include "spintrace_coul.h"
#include "boson_map.h"
#include <mathlink.h>


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
/// \brief Apply symmetry to Coulomb integrals (ab|cd) with a <= b, c <= d,
/// (a,b) <= (c,d) in reverse-lexicographical Bose map ordering (last entry more significant than first)
///
static void ApplySymmetry(boson_coords_t *ab, boson_coords_t *cd)
{
	// assume that all spatial orbitals are real,
	// so e.g. (ab|cd) = (ba|cd)
	if (ab[0] > ab[1]) { SwapBosonCoords(&ab[0], &ab[1]); }
	if (cd[0] > cd[1]) { SwapBosonCoords(&cd[0], &cd[1]); }

	// Coulomb operator is symmetric
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
/// \brief Generate the Coulomb integral (ab|cd) corresponding to spin-orbitals 'p0, p1, q0, q1';
/// output is encoded as bosonic bitfield
///
static inline void AddCoulIntegral(const fermi_coords_t p0, const fermi_coords_t p1, const fermi_coords_t q0, const fermi_coords_t q1)
{
	// spatial orbitals for Coulomb integrals (ab|cd)
	boson_coords_t ab[2] = { (boson_coords_t)(p0/2), (boson_coords_t)(q0/2) };	// round down
	boson_coords_t cd[2] = { (boson_coords_t)(p1/2), (boson_coords_t)(q1/2) };
	ApplySymmetry(ab, cd);

	// encode as bosons
	MLPutFunction(stdlink, "CoulInt", 2);
	// need ab[0] <= ab[1] and cd[0] <= cd[1]
	#ifndef NO_INT64_STRUCT
	MLPutInteger64(stdlink, (((bitfield_t)1) << ab[0]) + (((bitfield_t)1) << (ab[1]+1)));	// '+1' since we have bosons
	MLPutInteger64(stdlink, (((bitfield_t)1) << cd[0]) + (((bitfield_t)1) << (cd[1]+1)));
	#else
	// work-around for MLPutInteger64 on 32 bit systems
	bitfield_t w;
	w = (((bitfield_t)1) << ab[0]) + (((bitfield_t)1) << (ab[1]+1)); MLPutInteger64(stdlink, *(mlint64 *)&w);	// '+1' since we have bosons
	w = (((bitfield_t)1) << cd[0]) + (((bitfield_t)1) << (cd[1]+1)); MLPutInteger64(stdlink, *(mlint64 *)&w);
	#endif
}


//________________________________________________________________________________________________________________________
///
/// \brief Trace out spin and determine to-be calculated Coulomb integrals;
/// use fermionic coordinates as input, and encode spatial orbitals as bosonic bitfields
///
void SpinTraceCoulombBase(const fermi_coords_t p0, const fermi_coords_t p1, const fermi_coords_t q0, const fermi_coords_t q1)
{
	bool s0 = (p0-q0)%2 == 0 && (p1-q1)%2 == 0;	// whether both spins match
	bool s1 = (p0-q1)%2 == 0 && (p1-q0)%2 == 0;	// whether switched spins match

	if (s1)
	{
		if (s0)
		{
			MLPutFunction(stdlink, "Plus", 2);
			AddCoulIntegral(p0, p1, q0, q1);
		}
		// negative sign
		MLPutFunction(stdlink, "Times", 2);
		MLPutInteger(stdlink, -1);
		AddCoulIntegral(p0, p1, q1, q0);	// swap 'q0' and 'q1'
	}
	else	// s1 == false
	{
		if (s0)
		{
			AddCoulIntegral(p0, p1, q0, q1);
		}
		else
		{
			MLPutInteger(stdlink, 0);
		}
	}
}
