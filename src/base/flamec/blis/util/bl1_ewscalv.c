/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

void bl1_sewscalv( conj1_t conj, integer n, float* x, integer incx, float* y, integer incy )
{
	float*    chi;
	float*    psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_sscals( chi, psi );
	}
}

void bl1_dewscalv( conj1_t conj, integer n, double* x, integer incx, double* y, integer incy )
{
	double*   chi;
	double*   psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_dscals( chi, psi );
	}
}

void bl1_csewscalv( conj1_t conj, integer n, float* x, integer incx, scomplex* y, integer incy )
{
	float*    chi;
	scomplex* psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_csscals( chi, psi );
	}
}

void bl1_cewscalv( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy )
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

void bl1_zdewscalv( conj1_t conj, integer n, double* x, integer incx, dcomplex* y, integer incy )
{
	double*   chi;
	dcomplex* psi;
	integer       i;

	for ( i = 0; i < n; ++i )
	{
		chi = x + i*incx;
		psi = y + i*incy;

		bl1_zdscals( chi, psi );
	}
}

void bl1_zewscalv( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy )
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

