..  Copyright (C) 2024-2026, Advanced Micro Devices, Inc. All rights reserved.

..  Redistribution and use in source and binary forms, with or without
..  modification, are permitted provided that the following conditions are met:

..  1. Redistributions of source code must retain the above copyright notice,
..  this list of conditions and the following disclaimer.
..  2. Redistributions in binary form must reproduce the above copyright notice,
..  this list of conditions and the following disclaimer in the documentation
..  and/or other materials provided with the distribution.
..  3. Neither the name of the copyright holder nor the names of its
..  contributors may be used to endorse or promote products derived from this
..  software without specific prior written permission.

..  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
..  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
..  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
..  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
..  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
..  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
..  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
..  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
..  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
..  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
..  POSSIBILITY OF SUCH DAMAGE.

############
AOCL-LAPACK
############

AOCL-LAPACK is a high performant implementation of Linear Algebra PACKage (LAPACK).
LAPACK provides routines for solving systems of linear equations, least-squares problems,
eigenvalue problems, singular value problems, and the associated matrix factorizations. It is
extensible, easy to use, and available under an open-source license. Applications relying on standard
Netlib LAPACK interfaces can utilize AOCL-LAPACK with virtually no changes to their source
code. AOCL-LAPACK supports C, Fortran, and C++ template interfaces (for a subset of APIs) for
the LAPACK APIs.

AOCL-LAPACK is based on libFLAME, which was originally developed by current and former
members of the `Science of High-Performance Computing (SHPC) <https://shpc.oden.utexas.edu/>`_
group in the `Institute for Computational Engineering and Sciences <https://www.oden.utexas.edu/>`_
at `The University of Texas at Austin <https://www.utexas.edu/>`_ under the project name
libflame. The upstream libFLAME repository is available on `GitHub <https://github.com/flame/
libflame>`_. AMD is actively optimizing key routines in libFLAME as a part of the AOCL-LAPACK
library, for AMD “Zen”-based architectures in the "amd" fork of libFLAME hosted on AMD GitHub.

From AOCL 6.0, AOCL-LAPACK is conformant to LAPACK 3.12.1 specifications except for the Index-64
interfaces. In combination with the AOCL-BLAS library, which includes optimizations for
the AMD “Zen”-based processors, AOCL-LAPACK enables running high performing LAPACK
functionalities on AMD platforms.

AOCL-LAPACK supports C, Fortran, and C++ template interfaces for the LAPACK APIs.

.. toctree::
   :maxdepth: 1

   fortranInterfaces.rst

AOCL-LAPACK provides standard Fortran LAPACK APIs.
Please refer to `Netlib documentation <https://www.netlib.org/lapack/>`_ to see the API definitions.

.. toctree::
   :maxdepth: 1

   lapackeInterfaces.rst

AOCL-LAPACK provides standard C Interface to LAPACK APIs.
Please refer to `LAPACKE Netlib documentation <https://github.com/Reference-LAPACK/lapack/tree/master/LAPACKE>`_
to see the API definitions.

.. toctree::
   :maxdepth: 1

   cppInterfaces.rst

AOCL-LAPACK provides standard C++ Interface to LAPACK APIs. Please refer :ref:`CPP Interface documentation <cppInterfaces>`.

|

*AOCL-LAPACK is developed and maintained by AMD. For support, send an email to* toolchainsupport@amd.com.
