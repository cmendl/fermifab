/// \file spintrace_coul.h
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

#ifndef SPINTRACE_COUL_H
#define SPINTRACE_COUL_H

#include "fermi_map.h"


// trace out spin and determine to-be calculated Coulomb integrals
void SpinTraceCoulombBase(const fermi_coords_t p0, const fermi_coords_t p1, const fermi_coords_t q0, const fermi_coords_t q1);



#endif
