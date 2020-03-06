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

#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the sum of a list of integers
///
static inline int IntegerSum(const int *x, const int n)
{
	int sum = 0;
	int i;
	for (i = 0; i < n; i++)
	{
		sum += x[i];
	}

	return sum;
}


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the product of a list of integer numbers
///
static inline int IntegerProduct(const int *x, const int n)
{
	assert(n >= 0); // n == 0 is still reasonable

	int prod = 1;
	int i;
	for (i = 0; i < n; i++)
	{
		prod *= x[i];
	}

	return prod;
}


//________________________________________________________________________________________________________________________
//


// return 0 if k > n and 1 if k == 0
int Binomial(const int n, const int k);
