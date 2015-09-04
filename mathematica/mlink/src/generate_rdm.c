/// \file generate_rdm.c
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices.
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

#include "generate_rdm.h"
#include "fermi_map.h"
#include <mathlink.h>
#include <malloc.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices
/// on a Hilbert space with dimension 'orbs'
///
/// Given a normalized N-body wavefunction psi (coefficients of canonical Slater basis in lexicographical ordering),
/// <psi | K{i,j} psi> is the coefficient i,j of the p-body reduced density matrix of psi.
/// In the creation/annihilation formalism, K{i,j} is the operator a_j^dagger a_i acting on the N-body space wedge^N H.
///
void GenerateRDM(const int *orbs, const int *p1, const int *N1, const int *N2, const int nc)
{
	int i;

	const int N1tot = IntegerSum(N1, nc);
	const int N2tot = IntegerSum(N2, nc);
	const int p1tot = IntegerSum(p1, nc);
	const int p2tot = N1tot - N2tot + p1tot;

	// basic configuration setup
	fermi_config_t config;
	config.orbs = (int *)orbs;
	config.nc = nc;

	// fermionic particle spaces with configurations N1 and N2
	fermi_map_t baseMapN1, baseMapN2;
	config.N = (int *)N1; FermiMap(&config, &baseMapN1);
	config.N = (int *)N2; FermiMap(&config, &baseMapN2);

	// configurations p1 and p2
	fermi_map_t baseMapP1, baseMapP2;
	config.N = (int *)p1; FermiMap(&config, &baseMapP1);
	// N2 - p1 == N1 - p2
	int *p2 = (int *)malloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p2[i] = N1[i] - N2[i] + p1[i];
		assert(p2[i] >= 0);
	}
	config.N = p2; FermiMap(&config, &baseMapP2);

	// create sparse array
	MLPutFunction(stdlink, "SparseArray", 2);
	int numEntries = baseMapN2.num;
	for (i = 0; i < nc; i++)
	{
		numEntries *= Binomial(N2[i], p1[i])*Binomial(orbs[i]-N2[i]+p1[i], p2[i]);
	}
	MLPutFunction(stdlink, "List", numEntries);

	fermi_coords_t *y = (fermi_coords_t *)malloc(N1tot*sizeof(fermi_coords_t));

	int n[4];
	for (n[0] = 0; n[0] < baseMapP1.num; n[0]++)
	{
		for (n[1] = 0; n[1] < baseMapP2.num; n[1]++)
		{
			// fill first 'p2tot' coordinates of 'y' only
			FermiDecode(baseMapP2.map[n[1]], y, p2tot);

			for (n[3] = 0; n[3] < baseMapN2.num; n[3]++)
			{
				int sign[2];

				sign[0] = AnnihilSign(baseMapN2.map[n[3]], baseMapP1.map[n[0]]);
				if (!sign[0]) {
					continue;
				}
				// store remaining coordinates after annihilation in tail of 'y'
				FermiDecode(baseMapN2.map[n[3]] - baseMapP1.map[n[0]], y + p2tot, N2tot - p1tot);

				n[2] = Fermi2BaseSign(&baseMapN1, y, N1tot, &sign[1]);
				if (n[2] == -1) {
					continue;
				}

				// set entry
				MLPutFunction(stdlink, "Rule", 2);
				// Mathematica indices start at 1
				n[0]++; n[1]++; n[2]++; n[3]++;
				MLPutIntegerList(stdlink, n, 4);
				MLPutInteger(stdlink, sign[0]*sign[1]);
				// restore C++ indices
				n[0]--; n[1]--; n[2]--; n[3]--;
				numEntries--;	// keep track of number of entries for checking
			}
		}
	}
	assert(numEntries == 0);

	// set array dimensions
	n[0] = baseMapP1.num; n[1] = baseMapP2.num;
	n[2] = baseMapN1.num; n[3] = baseMapN2.num;
	MLPutIntegerList(stdlink, n, 4);

	// clean up
	free(y);
	free(p2);
	free(baseMapP2.map);
	free(baseMapP1.map);
	free(baseMapN2.map);
	free(baseMapN1.map);
}
