/// \file bitfield.h
/// \brief Bitfield representing fermionic occupations.
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

#ifndef BITFIELD_H
#define BITFIELD_H

#include "util.h"
#include <stdint.h>


//________________________________________________________________________________________________________________________
///
/// \brief Bitfield representing fermionic or bosonic occupations
///
typedef uint64_t bitfield_t;


int CompareBitfield(const void *x, const void *y);


//________________________________________________________________________________________________________________________
///
/// \brief Count the number of 1-bits in 'x'
///
///	Equivalent to the 'POPCNT' instruction (population count)
///
static inline int BitCount(bitfield_t x)
{
	int count = 0;

	while (x)
	{
		++count;
		x &= x - 1;	// remove last bit; also see the 'BLSR' instruction for calculating x & (x - 1)
	}

	return count;
}


//________________________________________________________________________________________________________________________
///
/// \brief Select the last trailing 1-bit
///
///	Equivalent to the 'BLSI' instruction (extract the lowest set isolated bit)
///
static inline bitfield_t LastBit(const bitfield_t x)
{
	return x & -x;	// -x: two's complement of x, e.g. x = 1110000_2 -> -x = ...1110010000_2
}


//________________________________________________________________________________________________________________________
///
/// \brief Sets trailing zeros to 1
///
///	Equivalent to the 'BLSFILL' instruction (fill from lowest set bit)
///
static inline bitfield_t BitTrailFill(const bitfield_t x)
{
	return x | (x - 1);		// e.g., ...011100_2 -> ...011111_2
}



#endif
