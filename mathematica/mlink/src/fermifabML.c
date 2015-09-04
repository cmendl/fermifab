/// \file fermifabML.c
/// \brief Main MathLink interface file.
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
#include "generate_rdm.h"
#include "spintrace_coul.h"
#include <mathlink.h>


//________________________________________________________________________________________________________________________
///
/// \brief MLSlaterRDM[s1_Integer, s2_Integer, p1_Integer]
///
void MLSlaterRDM(bitfield_t s1, bitfield_t s2, int p1)
{
	if (p1 <= 0)
	{
		MLEvaluate(stdlink, "Message[MLSlaterRDM::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	SlaterRDM(s1, s2, p1);
}


//________________________________________________________________________________________________________________________
///
/// \brief MLGenerateRDM[orbs_Integer, p_Integer, N_Integer]
///
void MLGenerateRDMs(int orbs, int p, int N)
{
	GenerateRDM(&orbs, &p, &N, &N, 1);
}

//________________________________________________________________________________________________________________________
///
/// \brief MLGenerateRDM[orbs_Integer, p1_Integer, N1_Integer, N2_Integer]
///
void MLGenerateRDMi(int orbs, int p1, int N1, int N2)
{
	GenerateRDM(&orbs, &p1, &N1, &N2, 1);
}

//________________________________________________________________________________________________________________________
///
/// \brief MLGenerateRDM[orbs_List, p1_List, N1_List, N2_List]
///
void MLGenerateRDMconf(int *orbs, long orbslen, int *p1, long p1len, int *N1, long N1len, int *N2, long N2len)
{
	if (orbslen != p1len || orbslen != N1len || orbslen != N2len)
	{
		MLEvaluate(stdlink, "Message[MLGenerateRDM::IncompatibleLength]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	GenerateRDM(orbs, p1, N1, N2, orbslen);
}


//________________________________________________________________________________________________________________________
///
/// \brief MLSpinTraceCoulombBase[p0_Integer, p1_Integer, q0_Integer, q1_Integer]
///
void MLSpinTraceCoulombBase(int p0, int p1, int q0, int q1)
{
	if (p0 == p1 || q0 == q1)
	{
		MLEvaluate(stdlink, "Message[MLSpinTraceCoulombBase::InvalidArg]");
		// discard 'ReturnPacket'
		MLNextPacket(stdlink);
		MLNewPacket(stdlink);	// discard
		// final output
		MLPutSymbol(stdlink, "$Failed");
		return;
	}

	// convert to zero-based indices
	SpinTraceCoulombBase(p0-1, p1-1, q0-1, q1-1);
}



#if MACINTOSH_MATHLINK

int main( int argc, char* argv[])
{
	/* Due to a bug in some standard C libraries that have shipped with
	 * MPW, zero is passed to MLMain below.  (If you build this program
	 * as an MPW tool, you can change the zero to argc.)
	 */
	argc = argc; /* suppress warning */
	return MLMain( 0, argv);
}

#elif WINDOWS_MATHLINK

#if __BORLANDC__
#pragma argsused
#endif

int PASCAL WinMain( HINSTANCE hinstCurrent, HINSTANCE hinstPrevious, LPSTR lpszCmdLine, int nCmdShow)
{
	char  buff[512];
	char FAR * buff_start = buff;
	char FAR * argv[32];
	char FAR * FAR * argv_end = argv + 32;

	hinstPrevious = hinstPrevious; /* suppress warning */

	if( !MLInitializeIcon( hinstCurrent, nCmdShow)) return 1;
	MLScanString( argv, &argv_end, &lpszCmdLine, &buff_start);
	return MLMain( (int)(argv_end - argv), argv);
}

#else

int main(int argc, char* argv[])
{
	return MLMain(argc, argv);
}

#endif
