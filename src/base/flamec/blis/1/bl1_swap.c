/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sswap( integer n, float* x, integer incx, float* y, integer incy )
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

void bl1_dswap( integer n, double* x, integer incx, double* y, integer incy )
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

void bl1_cswap( integer n, scomplex* x, integer incx, scomplex* y, integer incy )
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

void bl1_zswap( integer n, dcomplex* x, integer incx, dcomplex* y, integer incy )
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

