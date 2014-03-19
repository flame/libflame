/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sscal( int n, float* alpha, float* x, int incx )
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

void bl1_dscal( int n, double* alpha, double* x, int incx )
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

void bl1_csscal( int n, float* alpha, scomplex* x, int incx )
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

void bl1_cscal( int n, scomplex* alpha, scomplex* x, int incx )
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

void bl1_zdscal( int n, double* alpha, dcomplex* x, int incx )
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

void bl1_zscal( int n, dcomplex* alpha, dcomplex* x, int incx )
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

