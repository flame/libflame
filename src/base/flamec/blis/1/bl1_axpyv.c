/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_saxpyv( conj1_t conj, int n, float* alpha, float* x, int incx, float* y, int incy )
{
	bl1_saxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bl1_daxpyv( conj1_t conj, int n, double* alpha, double* x, int incx, double* y, int incy )
{
	bl1_daxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bl1_caxpyv( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy )
{
	scomplex* x_copy;
	int       incx_copy;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	x_copy    = x;
	incx_copy = incx;
	
	if ( bl1_is_conj( conj ) )
	{
		x_copy    = bl1_callocv( n );
		incx_copy = 1;
	
		bl1_ccopyv( conj,
		            n,
		            x,      incx,
		            x_copy, incx_copy );
	}

	bl1_caxpy( n,
	           alpha,
	           x_copy, incx_copy,
	           y,      incy );

	if ( bl1_is_conj( conj ) )
		bl1_cfree( x_copy );
}

void bl1_zaxpyv( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy )
{
	dcomplex* x_copy;
	int       incx_copy;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	x_copy    = x;
	incx_copy = incx;
	
	if ( bl1_is_conj( conj ) )
	{
		x_copy    = bl1_zallocv( n );
		incx_copy = 1;
	
		bl1_zcopyv( conj,
		            n,
		            x,      incx,
		            x_copy, incx_copy );
	}

	bl1_zaxpy( n,
	           alpha,
	           x_copy, incx_copy,
	           y,      incy );

	if ( bl1_is_conj( conj ) )
		bl1_zfree( x_copy );
}

