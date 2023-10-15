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

void bl1_snrm2( integer n, float* x, integer incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_snrm2( n,
	                     x, incx );
#else
	*norm = F77_snrm2( &n,
	                   x, &incx );
#endif
}

void bl1_dnrm2( integer n, double* x, integer incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dnrm2( n,
	                     x, incx );
#else
	*norm = F77_dnrm2( &n,
	                   x, &incx );
#endif
}

void bl1_cnrm2( integer n, scomplex* x, integer incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_scnrm2( n,
	                      x, incx );
#else
	*norm = F77_scnrm2( &n,
	                    x, &incx );
#endif
}

void bl1_znrm2( integer n, dcomplex* x, integer incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dznrm2( n,
	                      x, incx );
#else
	*norm = F77_dznrm2( &n,
	                    x, &incx );
#endif
}

