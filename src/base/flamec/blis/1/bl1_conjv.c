
#include "blis1.h"

void bl1_sconjv( int m, float* x, int incx )
{
	return;
}

void bl1_dconjv( int m, double* x, int incx )
{
	return;
}

void bl1_cconjv( int m, scomplex* x, int incx )
{
	float  m1        = bl1_sm1();
	float* x_conj    = ( float* ) x + 1;
	int    incx_conj = 2 * incx;

	bl1_sscal( m,
	           &m1,
	           x_conj, incx_conj );
}

void bl1_zconjv( int m, dcomplex* x, int incx )
{
	double  m1        = bl1_dm1();
	double* x_conj    = ( double* ) x + 1;
	int     incx_conj = 2 * incx;

	bl1_dscal( m,
	           &m1,
	           x_conj, incx_conj );
}

