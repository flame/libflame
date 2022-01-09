/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sewinvscalv( conj1_t conj, integer n, float* x, integer incx, float* y, integer incy )
{
	float*    chi;
	float*    psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_sinvscals( chi, psi );
	}
}

void bl1_dewinvscalv( conj1_t conj, integer n, double* x, integer incx, double* y, integer incy )
{
	double*   chi;
	double*   psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_dinvscals( chi, psi );
	}
}

void bl1_csewinvscalv( conj1_t conj, integer n, float* x, integer incx, scomplex* y, integer incy )
{
	float*    chi;
	scomplex* psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_csinvscals( chi, psi );
	}
}

void bl1_cewinvscalv( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy )
{
	scomplex* chi;
	scomplex* psi;
	scomplex  conjchi;
	integer       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_ccopyconj( chi, &conjchi );
			bl1_cinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_cinvscals( chi, psi );
		}
	}
}

void bl1_zdewinvscalv( conj1_t conj, integer n, double* x, integer incx, dcomplex* y, integer incy )
{
	double*   chi;
	dcomplex* psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_zdinvscals( chi, psi );
	}
}

void bl1_zewinvscalv( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy )
{
	dcomplex* chi;
	dcomplex* psi;
	dcomplex  conjchi;
	integer       i;

	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;

			bl1_zcopyconj( chi, &conjchi );
			bl1_zinvscals( &conjchi, psi );
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			chi = x + i*incx;
			psi = y + i*incy;
	
			bl1_zinvscals( chi, psi );
		}
	}
}

