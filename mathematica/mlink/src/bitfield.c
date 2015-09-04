/// \file bitfield.c
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

#include "bitfield.h"


//________________________________________________________________________________________________________________________
///
/// \brief Comparison function used for sorting
///
int CompareBitfield(const void *x, const void *y)
{
	const bitfield_t a = *(bitfield_t *)x;
	const bitfield_t b = *(bitfield_t *)y;

	if (a < b)
	{
		return -1;
	}
	else if (a == b)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}
