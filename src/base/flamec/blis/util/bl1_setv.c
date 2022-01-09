/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_isetv( integer n, integer* sigma, integer* x, integer incx )
{
	integer*   chi;
	integer    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_ssetv( integer n, float* sigma, float* x, integer incx )
{
	float* chi;
	integer    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_dsetv( integer n, double* sigma, double* x, integer incx )
{
	double* chi;
	integer     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		*chi = *sigma;
	}
}

void bl1_csetv( integer n, scomplex* sigma, scomplex* x, integer incx )
{
	scomplex* chi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

void bl1_zsetv( integer n, dcomplex* sigma, dcomplex* x, integer incx )
{
	dcomplex* chi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		chi->real = sigma->real;
		chi->imag = sigma->imag;
	}
}

