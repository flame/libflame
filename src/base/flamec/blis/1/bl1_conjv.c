/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sconjv( integer m, float* x, integer incx )
{
	return;
}

void bl1_dconjv( integer m, double* x, integer incx )
{
	return;
}

void bl1_cconjv( integer m, scomplex* x, integer incx )
{
	float  m1        = bl1_sm1();
	float* x_conj    = ( float* ) x + 1;
	integer    incx_conj = 2 * incx;

	bl1_sscal( m,
	           &m1,
	           x_conj, incx_conj );
}

void bl1_zconjv( integer m, dcomplex* x, integer incx )
{
	double  m1        = bl1_dm1();
	double* x_conj    = ( double* ) x + 1;
	integer     incx_conj = 2 * incx;

	bl1_dscal( m,
	           &m1,
	           x_conj, incx_conj );
}

