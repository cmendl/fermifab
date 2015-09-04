FermiFab toolbox for fermionic many-particle quantum systems
============================================================

The FermiFab toolbox consists of a Matlab and Mathematica implementation, which are independent of each other but use the same principles and basically perform the same calculations. Both are supported by C interfaces (MEX and MathLink, respectively) to speed up computations and provide an additional abstraction layer.

Matlab
------
- Installation: add the matlab/fermifab subfolder to your Matlab search path
- MEX binaries: if your platform-specific MEX binaries in matlab/fermifab happen to be missing: make sure MEX has been properly configured for your system (`mex -setup`), call `build_mex.m` and hope that everything compiles well!
- Folder structure:
    * fermifab:         Matlab toolbox files
    * mex:              MEX C source code
    * doc:              documentation
    * examples:         example and demonstration files
    * test:             testing and verification

Mathematica
-----------
- Demo file illustrating the reduced-density matrix (RDM) formalism: `RDMdemo.nb`
- The main toolbox code resides in `fermifab.m` and the MathLink C backend
- The 'mlink/fermifabML' subdirectory contains binary MathLink executables for several platforms; if these are not compatible with your architecture, consult the local readme file for details how to compile the MathLink source code

License
-------
Copyright (c) 2008-2015, Christian B. Mendl  
All rights reserved.  
http://christian.mendl.net

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License  
http://www.opensource.org/licenses/bsd-license.php

Reference
---------
  Christian B. Mendl  
  The FermiFab toolbox for fermionic many-particle quantum systems  
  Comput. Phys. Commun. 182, 1327-1337 (2011) [DOI](http://dx.doi.org/10.1016/j.cpc.2011.01.028)  
  preprint http://arxiv.org/abs/1103.0872
