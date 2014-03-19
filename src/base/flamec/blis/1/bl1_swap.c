/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sswap( int n, float* x, int incx, float* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_sswap( n,
	             x, incx, 
	             y, incy );
#else
	F77_sswap( &n,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_dswap( int n, double* x, int incx, double* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_dswap( n,
	             x, incx, 
	             y, incy );
#else
	F77_dswap( &n,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_cswap( int n, scomplex* x, int incx, scomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_cswap( n,
	             x, incx, 
	             y, incy );
#else
	F77_cswap( &n,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_zswap( int n, dcomplex* x, int incx, dcomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zswap( n,
	             x, incx, 
	             y, incy );
#else
	F77_zswap( &n,
	           x, &incx, 
	           y, &incy );
#endif
}

