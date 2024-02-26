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

void bl1_saxpyv( conj1_t conj, integer n, float* alpha, float* x, integer incx, float* y, integer incy )
{
	bl1_saxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bl1_daxpyv( conj1_t conj, integer n, double* alpha, double* x, integer incx, double* y, integer incy )
{
	bl1_daxpy( n,
	           alpha,
	           x, incx,
	           y, incy );
}

void bl1_caxpyv( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy )
{
	scomplex* x_copy;
	integer       incx_copy;

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

void bl1_zaxpyv( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy )
{
	dcomplex* x_copy;
	integer       incx_copy;

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

