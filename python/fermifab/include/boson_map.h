/// \file boson_map.h
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

#pragma once

#include "bitfield.h"


//________________________________________________________________________________________________________________________
///
/// \brief Bosonic 'coordinates', i.e., orbital (or mode) of each particle 
///
typedef int boson_coords_t;


// encode coordinates 'x' into a bitfield
bitfield_t BosonEncode(const boson_coords_t *x, const int N);

// decode the bitfield 'w' into coordinates 'x'
void BosonDecode(const bitfield_t w, boson_coords_t *x, const int N);


//________________________________________________________________________________________________________________________
///
/// \brief Bose map for enumeration of bit-encoded occupations
///
typedef struct
{
    bitfield_t *map;	//!< list of bit-encoded occupations
    int num;			//!< length of the list
}
boson_map_t;


// map base indices to bit-encoded coordinates
int BosonMap(const int modes, const int N, boson_map_t *bm);

// convert bosonic coordinates to base index
int Boson2Base(const boson_map_t *bm, const boson_coords_t *x, const int N);
