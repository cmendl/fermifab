/// \file fermi_map.h
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

#ifndef FERMI_MAP_H
#define FERMI_MAP_H

#include "bitfield.h"


//________________________________________________________________________________________________________________________
///
/// \brief Fermionic 'coordinates', i.e., list of occupied orbitals
///
typedef int fermi_coords_t;


// encode coordinates 'x' into a bitfield
bitfield_t FermiEncode(const fermi_coords_t *x, const int N);

// decode the bitfield 'w' into coordinates 'x'
void FermiDecode(const bitfield_t w, fermi_coords_t *x, const int N);


//________________________________________________________________________________________________________________________
///
/// \brief Fermi Hilbert space configuration, i.e., number of orbitals and particles
///
typedef struct
{
	int *orbs;		//!< number of orbitals in each partition
	int *N;			//!< corresponding number of particles
	int nc;			//!< number of partitions
}
fermi_config_t;


//________________________________________________________________________________________________________________________
///
/// \brief Fermi map for enumeration of bit-encoded occupations
///
typedef struct
{
	bitfield_t *map;	//!< list of bit-encoded occupations
	int num;			//!< length of the list
}
fermi_map_t;


// map base indices to bit-encoded coordinates
void FermiMap(const fermi_config_t *config, fermi_map_t *fm);

// lexicographically next bit permutation
bitfield_t NextFermi(const bitfield_t f);

// convert Fermi coordinates to base index and permutation sign
int Fermi2Base    (const fermi_map_t *fm, const fermi_coords_t *x, const int N);
int Fermi2BaseSign(const fermi_map_t *fm, const fermi_coords_t *x, const int N, int *sign);


// annihilation sign
int AnnihilSign(const bitfield_t n, const bitfield_t a);

// annihilation sign mask
bitfield_t AnnihilSignMask(const bitfield_t f);


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the permutation sign of (length,length-1,...,1)
///
static inline int ReversePermSign(const int length)
{
	int n = (length-1)*length/2;
	return 1 - 2*(n % 2);
}



#endif
