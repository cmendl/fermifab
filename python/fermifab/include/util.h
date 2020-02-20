/// \file util.h
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

#pragma once


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the sum of a list of integers
///
static inline int IntegerSum(const int *list, const int n)
{
	int i, ret;

	ret = 0;
	for (i = 0; i < n; i++)
	{
		ret += list[i];
	}

	return ret;
}


// return 0 if k > n and 1 if k == 0
int Binomial(const int n, const int k);
