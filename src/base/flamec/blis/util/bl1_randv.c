/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_srandv( integer n, float* x, integer incx )
{
	float* chi;
	integer    i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_srands( chi );
	}
}

void bl1_drandv( integer n, double* x, integer incx )
{
	double* chi;
	integer     i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_drands( chi );
	}
}

void bl1_crandv( integer n, scomplex* x, integer incx )
{
	scomplex* chi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_crands( chi );
	}
}

void bl1_zrandv( integer n, dcomplex* x, integer incx )
{
	dcomplex* chi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;

		bl1_zrands( chi );
	}
}

