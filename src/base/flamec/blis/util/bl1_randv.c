/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

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

