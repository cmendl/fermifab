FermiFab toolbox for fermionic many-particle quantum systems
============================================================

FermiFab is a quantum physics toolbox for small fermionic systems. It can represent fermionic many-body quantum states and calculate their reduced density matrices (RDMs). The toolbox transparently handles the inherent antisymmetrization of wavefunctions and incorporates the creation and annihilation formalism. FermiFab offers independent implementations in Python, Matlab and Mathematica, which are supported by C interfaces to speed up computations and provide an additional abstraction layer.

Python
------
![](https://github.com/cmendl/fermifab/workflows/python_unittest/badge.svg)

- Installation:
  1. requirements: `sudo apt install libblas-dev liblapacke-dev` (on Ubuntu Linux) or similar, and `pip install -r python/requirements.txt`
  2. run `python setup.py install --user` from the python/ subfolder (or `python3 ...` instead of `python ...`)
- Documentation: [fermifab.readthedocs.io](https://fermifab.readthedocs.io), [user guide](python/doc/user_guide.ipynb)

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
Copyright (c) 2008-2020, Ismael Medina Suárez and Christian B. Mendl  
All rights reserved.  

This program is free software; you can redistribute it and/or
modify it under the terms of the Simplified BSD License  
http://www.opensource.org/licenses/bsd-license.php

Reference
---------
  Christian B. Mendl  
  The FermiFab toolbox for fermionic many-particle quantum systems  
  Comput. Phys. Commun. 182, 1327-1337 (2011) [DOI](http://dx.doi.org/10.1016/j.cpc.2011.01.028)  
  preprint http://arxiv.org/abs/1103.0872
