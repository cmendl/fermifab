/// \file util.h
/// \brief Utility functions.
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

#ifndef UTIL_H
#define UTIL_H

#include <stdbool.h>

#if defined(_MSC_VER) && (_MSC_VER <= 1600)

// 'bool' not known to older versions of Visual Studio C
typedef int bool;

#define true	1
#define false	0

#endif


// Microsoft specific 'inline'
#if defined(_MSC_VER)
#define inline __inline
#endif


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the minimum of two integers
///
static inline int mini(const int x, const int y)
{
	if (x < y)
	{
		return x;
	}
	else
	{
		return y;
	}
}

//________________________________________________________________________________________________________________________
///
/// \brief Calculate the maximum of two integers
///
static inline int maxi(const int x, const int y)
{
	if (y < x)
	{
		return x;
	}
	else
	{
		return y;
	}
}


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



#endif
