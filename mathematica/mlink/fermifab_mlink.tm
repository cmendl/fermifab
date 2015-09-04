
:Begin:
:Function:		MLSlaterRDM
:Pattern:		MLSlaterRDM[s1_Integer, s2_Integer, p1_Integer]
:Arguments:		{ s1, s2, p1 }
:ArgumentTypes:	{ Integer64, Integer64, Integer }
:ReturnType:	Manual
:End:

:Evaluate: MLSlaterRDM::usage = "MLSlaterRDM[s1_Integer,s2_Integer,p1_Integer] calculates the p1-particle reduced density matrices (RDMs) from Fermi states (Slater determinants) s1 and s2."
:Evaluate: MLSlaterRDM::InvalidArg = "Syntax: MLSlaterRDM[s1_Integer,s2_Integer,p1_Integer], p1 must be positive."



:Begin:
:Function:		MLGenerateRDMs
:Pattern:		MLGenerateRDM[orbs_Integer, p_Integer, N_Integer]
:Arguments:		{ orbs, p, N }
:ArgumentTypes:	{ Integer, Integer, Integer }
:ReturnType:	Manual
:End:

:Begin:
:Function:		MLGenerateRDMi
:Pattern:		MLGenerateRDM[orbs_Integer, p1_Integer, N1_Integer, N2_Integer]
:Arguments:		{ orbs, p1, N1, N2 }
:ArgumentTypes:	{ Integer, Integer, Integer, Integer }
:ReturnType:	Manual
:End:

:Begin:
:Function:		MLGenerateRDMconf
:Pattern:		MLGenerateRDM[orbs_List, p1_List, N1_List, N2_List]
:Arguments:		{ orbs, p1, N1, N2 }
:ArgumentTypes:	{ IntegerList, IntegerList, IntegerList, IntegerList }
:ReturnType:	Manual
:End:

:Evaluate: MLGenerateRDM::usage = "MLGenerateRDM[orbs,p1,N1,N2] generates the RDM kernel for p-body reduced density matrices."
:Evaluate: MLGenerateRDM::IncompatibleLength = "All input arguments must be of the same length."



:Begin:
:Function:		MLSpinTraceCoulombBase
:Pattern:		MLSpinTraceCoulombBase[p0_Integer, p1_Integer, q0_Integer, q1_Integer]
:Arguments:		{ p0, p1, q0, q1 }
:ArgumentTypes:	{ Integer, Integer, Integer, Integer }
:ReturnType:	Manual
:End:

:Evaluate: MLSpinTraceCoulombBase::usage = "MLSpinTraceCoulombBase[p0_Integer, p1_Integer, q0_Integer, q1_Integer] traces out spin and determines to-be calculated Coulomb integrals."
:Evaluate: MLSpinTraceCoulombBase::InvalidArg = "Syntax: MLSpinTraceCoulombBase[p0_Integer, p1_Integer, q0_Integer, q1_Integer], (p0,p1) and (q0,q1) must be pairwise different."
