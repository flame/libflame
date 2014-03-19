/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_icopyv( conj1_t conj, int m, int* x, int incx, int* y, int incy )
{
	int*      chi;
	int*      psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

void bl1_scopyv( conj1_t conj, int m, float* x, int incx, float* y, int incy )
{
	bl1_scopy( m,
	           x, incx, 
	           y, incy );
}

void bl1_dcopyv( conj1_t conj, int m, double* x, int incx, double* y, int incy )
{
	bl1_dcopy( m,
	           x, incx, 
	           y, incy );
}

void bl1_ccopyv( conj1_t conj, int m, scomplex* x, int incx, scomplex* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	bl1_ccopy( m,
	           x, incx, 
	           y, incy );

	if ( bl1_is_conj( conj ) )
		bl1_cconjv( m,
	                y, incy );
}

void bl1_zcopyv( conj1_t conj, int m, dcomplex* x, int incx, dcomplex* y, int incy )
{
	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	bl1_zcopy( m,
	           x, incx, 
	           y, incy );

	if ( bl1_is_conj( conj ) )
		bl1_zconjv( m,
	                y, incy );
}

// --- Mixed-datatype and general stride copy routines---------------

// sd ds
void bl1_sdcopyv( conj1_t conj, int m, float* x, int incx, double* y, int incy )
{
	float*    chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}
void bl1_dscopyv( conj1_t conj, int m, double* x, int incx, float* y, int incy )
{
	double*   chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = *chi;

		chi += incx;
		psi += incy;
	}
}

// sc cs
void bl1_sccopyv( conj1_t conj, int m, float* x, int incx, scomplex* y, int incy )
{
	float*    chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bl1_cscopyv( conj1_t conj, int m, scomplex* x, int incx, float* y, int incy )
{
	scomplex* chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// sz zs
void bl1_szcopyv( conj1_t conj, int m, float* x, int incx, dcomplex* y, int incy )
{
	float*    chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bl1_zscopyv( conj1_t conj, int m, dcomplex* x, int incx, float* y, int incy )
{
	dcomplex* chi;
	float*    psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dc cd
void bl1_dccopyv( conj1_t conj, int m, double* x, int incx, scomplex* y, int incy )
{
	double*   chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0F;

		chi += incx;
		psi += incy;
	}
}
void bl1_cdcopyv( conj1_t conj, int m, scomplex* x, int incx, double* y, int incy )
{
	scomplex* chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// dz zd
void bl1_dzcopyv( conj1_t conj, int m, double* x, int incx, dcomplex* y, int incy )
{
	double*   chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = *chi;
		psi->imag = 0.0;

		chi += incx;
		psi += incy;
	}
}
void bl1_zdcopyv( conj1_t conj, int m, dcomplex* x, int incx, double* y, int incy )
{
	dcomplex* chi;
	double*   psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		*psi = chi->real;

		chi += incx;
		psi += incy;
	}
}

// cz zc
void bl1_czcopyv( conj1_t conj, int m, scomplex* x, int incx, dcomplex* y, int incy )
{
	scomplex* chi;
	dcomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bl1_is_conj( conj ) )
		bl1_zconjv( m,
	                y, incy );
}
void bl1_zccopyv( conj1_t conj, int m, dcomplex* x, int incx, scomplex* y, int incy )
{
	dcomplex* chi;
	scomplex* psi;
	int       i;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// Initialize pointers.
	chi = x;
	psi = y;

	for ( i = 0; i < m; ++i )
	{
		psi->real = chi->real;
		psi->imag = chi->imag;

		chi += incx;
		psi += incy;
	}

	if ( bl1_is_conj( conj ) )
		bl1_cconjv( m,
	                y, incy );
}

