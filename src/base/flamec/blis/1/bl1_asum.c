/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

void bl1_sasum( integer n, float* x, integer incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_sasum( n,
	                     x, incx );
#else
	*norm = F77_sasum( &n,
	                   x, &incx );
#endif
}

void bl1_dasum( integer n, double* x, integer incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dasum( n,
	                     x, incx );
#else
	*norm = F77_dasum( &n,
	                   x, &incx );
#endif
}

void bl1_casum( integer n, scomplex* x, integer incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_scasum( n,
	                      x, incx );
#else
	*norm = F77_scasum( &n,
	                    x, &incx );
#endif
}

void bl1_zasum( integer n, dcomplex* x, integer incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dzasum( n,
	                      x, incx );
#else
	*norm = F77_dzasum( &n,
	                    x, &incx );
#endif
}

