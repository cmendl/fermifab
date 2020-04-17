/// \file generate_rdm.c
/// \brief Create the matrix kernels K required for calculating p-body reduced density matrices.
//
//  Copyright (c) 2008-2020, Christian B. Mendl
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

#include "generate_rdm.h"
#include "fermi_map.h"
#include "util.h"
#include <malloc.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Generate the kernel tensor K required for calculating p-body reduced density matrices
/// on a Hilbert space with dimension 'orbs'
///
/// Given a normalized N-body wavefunction psi (coefficients of canonical Slater basis in lexicographical ordering),
/// <psi | K{i,j} psi> is the coefficient i,j of the p-body reduced density matrix of psi.
/// In the creation/annihilation formalism, K{i,j} is the operator a_j^dagger a_i acting on the N-body space wedge^N H.
///
int GenerateRDM(const int *orbs, const int *p1, const int *N1, const int *N2, const int nc, sparse_array_t *K)
{
	int i;
	int status;

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
	config.N = (int *)N1; status = FermiMap(&config, &baseMapN1); if (status < 0) { return status; }
	config.N = (int *)N2; status = FermiMap(&config, &baseMapN2); if (status < 0) { return status; }

	// configurations p1 and p2
	fermi_map_t baseMapP1, baseMapP2;
	config.N = (int *)p1; status = FermiMap(&config, &baseMapP1); if (status < 0) { return status; }
	// N2 - p1 == N1 - p2
	int *p2 = (int *)malloc(nc*sizeof(int));
	for (i = 0; i < nc; i++)
	{
		p2[i] = N1[i] - N2[i] + p1[i];
		assert(p2[i] >= 0);
	}
	config.N = p2; status = FermiMap(&config, &baseMapP2); if (status < 0) { return status; }

	// create sparse array
	K->rank = 4;
	K->dims = (int *)malloc(K->rank * sizeof(int));
	K->nnz = baseMapN2.num;
	for (i = 0; i < nc; i++)
	{
		K->nnz *= Binomial(N2[i], p1[i])*Binomial(orbs[i]-N2[i]+p1[i], p2[i]);
	}
	K->val = (double *)malloc(K->nnz * sizeof(double));    if (K->val == NULL) { return -1; }
	K->ind = (int *)malloc(K->nnz*K->rank * sizeof(int));  if (K->ind == NULL) { return -1; }

	fermi_coords_t *y = (fermi_coords_t *)malloc(N1tot*sizeof(fermi_coords_t));

	int count = 0;
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
				K->ind[4*count  ] = n[0];
				K->ind[4*count+1] = n[1];
				K->ind[4*count+2] = n[2];
				K->ind[4*count+3] = n[3];
				K->val[count] = sign[0]*sign[1];
				count++;
			}
		}
	}
	assert(count == K->nnz);

	// set array dimensions
	K->dims[0] = baseMapP1.num;
	K->dims[1] = baseMapP2.num;
	K->dims[2] = baseMapN1.num;
	K->dims[3] = baseMapN2.num;

	// clean up
	free(y);
	free(p2);
	free(baseMapP2.map);
	free(baseMapP1.map);
	free(baseMapN2.map);
	free(baseMapN1.map);

	return 0;
}
