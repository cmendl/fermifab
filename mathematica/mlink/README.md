C source code of the FermiFab toolbox (MathLink interface)
==========================================================

The C source code should compile on most platforms without modifications; the binary executables are then stored in the 'mlink/fermifabML' subdirectory.

To compile the code under
- Windows (using, e.g., Microsoft Visual Studio): set up a project file, add the C source code files in the 'include' and 'src' subdirectories, provide the path to `mathlink.h` from the Mathematica installation and add the required Mathematica library files as described in the Mathematica MathLink documentation
- Linux, Unix etc: `make` should work after adapting the paths in `Makefile`
- other platforms: adapt the makefile to your system; see also the Mathematica MathLink documentation
