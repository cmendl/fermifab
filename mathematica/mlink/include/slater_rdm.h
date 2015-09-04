/// \file slater_rdm.h
/// \brief Calculate the reduced density matrix of the outer product of two Slater determinants.
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

#ifndef SLATER_RDM_H
#define SLATER_RDM_H

#include "bitfield.h"


// calculate the reduced density matrix of the outer product of two Slater determinants
void SlaterRDM(const bitfield_t s1, const bitfield_t s2, const int p1);



#endif
