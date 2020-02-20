/// \file util.c
/// \brief Utility functions.
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

#include "util.h"
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the binomial coefficient 'n choose k' recursively
///
int Binomial(const int n, const int k)
{
	assert(n >= 0 && k >= 0);

	if (k > n)
	{
		return 0;
	}
	else if (k == 0)
	{
		return 1;
	}

	// recursive function call
	return Binomial(n-1, k) + Binomial(n-1, k-1);
}
