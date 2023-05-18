# AOCL-LAPACK library

AOCL-LAPACK is a high performant implementation of Linear Algebra PACKage(LAPACK).
LAPACK provides routines for solving systems of linear equations, least-squares 
problems, eigenvalue problems, singular value problems, and the associated 
matrix factorizations. It is extensible, easy to use, and available under 
open-source BSD-3 license. AOCL-LAPACK is a C-only implementation. Applications 
relying on standard Netlib LAPACK interfaces can utilize libFLAME with virtually 
no changes to their source code. 

AOCL-LAPACK was originally developed by current and former members of the 
[Science of High-Performance Computing](http://shpc.ices.utexas.edu/)
(SHPC) group in the
[Institute for Computational Engineering and Sciences](https://www.ices.utexas.edu/)
at [The University of Texas at Austin](https://www.utexas.edu/) under project 
name libflame. The upstream libflame repository is available 
[here](https://github.com/flame/libflame). AMD is actively optimizing key routines 
in libFLAME as part of AOCL-LAPACK library, for AMD Zen core based architectures 
in the "amd" fork of libFLAME hosted on [github](https://github.com/amd/libflame).

For detailed instructions on how to configure, build, install, and link against 
libflame on AMD CPUs, please refer to the AOCL User Guide located on AMD 
developer [portal](https://www.amd.com/en/developer/aocl.html).

Upstream repository contains libflame reference manual and a complete API 
reference. If you have LaTeX installed on your system, you may simply change 
into the 'docs/libflame' subdirectory of the top-level directory of the 
libflame source tree and build the document from its source. You may also 
find a copy of the document [here](docs/libflame/libflame.pdf) on github.

For any issues/suggestion in the "amd" fork of libFLAME, please email 
toolchainsupport@amd.com

Also, please read the LICENSE file for information on copying and distributing 
this software.

