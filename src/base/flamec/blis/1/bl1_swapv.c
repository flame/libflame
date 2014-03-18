
#include "blis1.h"

void bl1_sswapv( int n, float* x, int incx, float* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_sswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_dswapv( int n, double* x, int incx, double* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_dswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_cswapv( int n, scomplex* x, int incx, scomplex* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_cswap( n,
	           x, incx, 
	           y, incy );
}

void bl1_zswapv( int n, dcomplex* x, int incx, dcomplex* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	bl1_zswap( n,
	           x, incx, 
	           y, incy );
}

