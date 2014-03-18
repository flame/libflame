
#include "blis1.h"

void bl1_srandv( int n, float* x, int incx )
{
	float* chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_srands( chi );
	}
}

void bl1_drandv( int n, double* x, int incx )
{
	double* chi;
	int     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_drands( chi );
	}
}

void bl1_crandv( int n, scomplex* x, int incx )
{
	scomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_crands( chi );
	}
}

void bl1_zrandv( int n, dcomplex* x, int incx )
{
	dcomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_zrands( chi );
	}
}

