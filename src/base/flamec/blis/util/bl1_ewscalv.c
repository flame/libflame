/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sewscalv( conj1_t conj, int n, float* x, int incx, float* y, int incy )
{
	float*    chi;
	float*    psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_sscals( chi, psi );
	}
}

void bl1_dewscalv( conj1_t conj, int n, double* x, int incx, double* y, int incy )
{
	double*   chi;
	double*   psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_dscals( chi, psi );
	}
}

void bl1_csewscalv( conj1_t conj, int n, float* x, int incx, scomplex* y, int incy )
{
	float*    chi;
	scomplex* psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_csscals( chi, psi );
	}
}

void bl1_cewscalv( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy )
{
	scomplex* chi;
	scomplex* psi;
	scomplex  conjchi;
	int       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_ccopyconj( chi, &conjchi );
			bl1_cscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_cscals( chi, psi );
		}
	}
}

void bl1_zdewscalv( conj1_t conj, int n, double* x, int incx, dcomplex* y, int incy )
{
	double*   chi;
	dcomplex* psi;
	int       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_zdscals( chi, psi );
	}
}

void bl1_zewscalv( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy )
{
	dcomplex* chi;
	dcomplex* psi;
	dcomplex  conjchi;
	int       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_zcopyconj( chi, &conjchi );
			bl1_zscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_zscals( chi, psi );
		}
	}
}

