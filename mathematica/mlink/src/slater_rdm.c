/// \file slater_rdm.c
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

#include "slater_rdm.h"
#include "fermi_map.h"
#include <mathlink.h>
#include <assert.h>


//________________________________________________________________________________________________________________________
///
/// \brief Calculate the p-particle reduced density matrices (RDMs)
/// of the outer product of two fermionic Slater determinants 's1' and 's2': |s1Xs2|
///
/// Encoding: orbital 0 corresponds to LSB (least significant bit)
///
void SlaterRDM(const bitfield_t s1, const bitfield_t s2, const int p1)
{
	// number of particles in 's1'
	int n1 = BitCount(s1);
	if (p1 > n1)
	{
		MLPutFunction(stdlink, "List", 0);
	}
	else
	{
		// "force" mask
		bitfield_t fmask = s1 ^ s2;

		// states which must be annihilated in 's1'
		bitfield_t sforce1 = fmask & s1;

		// number of to-be annihilated "choice" states in 's1'
		int nchoice1 = p1 - BitCount(sforce1);
		if (nchoice1 < 0)
		{
			MLPutFunction(stdlink, "List", 0);
		}
		else
		{
			bitfield_t sforce2 = fmask & s2;

			// annihilation sign masks
			bitfield_t amask1 = AnnihilSignMask(s1) << 1;
			bitfield_t amask2 = AnnihilSignMask(s2) << 1;

			int factor = ReversePermSign(p1)*ReversePermSign(BitCount(s2)-n1+p1)*
				IntegerParitySign(amask1 & sforce1)*IntegerParitySign(amask2 & sforce2);

			if (nchoice1 == 0)
			{
				// no choice states

				// create list as return value
				MLPutFunction(stdlink, "List", 1);

				MLPutFunction(stdlink, "Times", 2);
				MLPutInteger(stdlink, factor);
				MLPutFunction(stdlink, "IXI", 2);
				#ifndef NO_INT64_STRUCT
				MLPutInteger64(stdlink, sforce1);
				MLPutInteger64(stdlink, sforce2);
				#else
				MLPutInteger64(stdlink, *(mlint64 *)&sforce1);
				MLPutInteger64(stdlink, *(mlint64 *)&sforce2);
				#endif
			}
			else
			{
				// nchoice1 > 0

				// "choice" mask
				bitfield_t cmask = s1 & s2;

				// number of "choice" orbitals
				int kchoice = BitCount(cmask);
				assert(nchoice1 <= kchoice);

				// create list as return value
				MLPutFunction(stdlink, "List", Binomial(kchoice, nchoice1));

				// iterate Fermi map for 'choice' orbs
				bitfield_t t = (((bitfield_t)1) << nchoice1) - 1;
				do
				{
					bitfield_t schoice = BitDistribute(t, cmask);

					MLPutFunction(stdlink, "Times", 2);
					MLPutInteger(stdlink, factor*IntegerParitySign(amask1 & schoice)*IntegerParitySign(amask2 & schoice));
					MLPutFunction(stdlink, "IXI", 2);
					#ifndef NO_INT64_STRUCT
					MLPutInteger64(stdlink, sforce1 + schoice);
					MLPutInteger64(stdlink, sforce2 + schoice);
					#else
					// work-around for MLPutInteger64 on 32 bit systems
					bitfield_t w;
					w = sforce1 + schoice; MLPutInteger64(stdlink, *(mlint64 *)&w);
					w = sforce2 + schoice; MLPutInteger64(stdlink, *(mlint64 *)&w);
					#endif

					// next Fermi state
					t = NextFermi(t);
				}
				while ((t >> kchoice) == 0);
			}
		}
	}
}
