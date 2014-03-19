/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_isetv( int n, int* sigma, int* x, int incx )
{
	int*   chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_ssetv( int n, float* sigma, float* x, int incx )
{
	float* chi;
	int    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_dsetv( int n, double* sigma, double* x, int incx )
{
	double* chi;
	int     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_csetv( int n, scomplex* sigma, scomplex* x, int incx )
{
	scomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

void bl1_zsetv( int n, dcomplex* sigma, dcomplex* x, int incx )
{
	dcomplex* chi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

