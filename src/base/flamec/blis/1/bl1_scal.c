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

void bl1_sscal( integer n, float* alpha, float* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_sscal( n,
	             *alpha,
	             x, incx );
#else
	F77_sscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_dscal( integer n, double* alpha, double* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_dscal( n,
	             *alpha,
	             x, incx );
#else
	F77_dscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_csscal( integer n, float* alpha, scomplex* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_csscal( n,
	              *alpha,
	              x, incx );
#else
	F77_csscal( &n,
	            alpha,
	            x, &incx );
#endif
}

void bl1_cscal( integer n, scomplex* alpha, scomplex* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_cscal( n,
	             alpha,
	             x, incx );
#else
	F77_cscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_zdscal( integer n, double* alpha, dcomplex* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zdscal( n,
	              *alpha,
	              x, incx );
#else
	F77_zdscal( &n,
	            alpha,
	            x, &incx );
#endif
}

void bl1_zscal( integer n, dcomplex* alpha, dcomplex* x, integer incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zscal( n,
	             alpha,
	             x, incx );
#else
	F77_zscal( &n,
	           alpha,
	           x, &incx );
#endif
}

