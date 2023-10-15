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

void bl1_sapdiagmv( side1_t side, conj1_t conj, integer m, integer n, float* x, integer incx, float* a, integer a_rs, integer a_cs )
{
	float*    chi;
	float*    a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_sewscalv( conj,
			              n_elem,
			              x,       incx,
			              a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_sscalv( conj,
			            n_elem,
			            chi,
			            a_begin, inca );
		}
	}
}

void bl1_dapdiagmv( side1_t side, conj1_t conj, integer m, integer n, double* x, integer incx, double* a, integer a_rs, integer a_cs )
{
	double*   chi;
	double*   a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_dewscalv( conj,
			              n_elem,
			              x,       incx,
			              a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_dscalv( conj,
			            n_elem,
			            chi,
			            a_begin, inca );
		}
	}
}

void bl1_csapdiagmv( side1_t side, conj1_t conj, integer m, integer n, float* x, integer incx, scomplex* a, integer a_rs, integer a_cs )
{
	float*    chi;
	scomplex* a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_csewscalv( conj,
			               n_elem,
			               x,       incx,
			               a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_csscalv( conj,
			             n_elem,
			             chi,
			             a_begin, inca );
		}
	}
}

void bl1_capdiagmv( side1_t side, conj1_t conj, integer m, integer n, scomplex* x, integer incx, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* chi;
	scomplex* a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_cewscalv( conj,
			              n_elem,
			              x,       incx,
			              a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_cscalv( conj,
			            n_elem,
			            chi,
			            a_begin, inca );
		}
	}
}

void bl1_zdapdiagmv( side1_t side, conj1_t conj, integer m, integer n, double* x, integer incx, dcomplex* a, integer a_rs, integer a_cs )
{
	double*   chi;
	dcomplex* a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_zdewscalv( conj,
			               n_elem,
			               x,       incx,
			               a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_zdscalv( conj,
			             n_elem,
			             chi,
			             a_begin, inca );
		}
	}
}

void bl1_zapdiagmv( side1_t side, conj1_t conj, integer m, integer n, dcomplex* x, integer incx, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* chi;
	dcomplex* a_begin;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	// An optimization: if A is row-major, then we can proceed as if the
	// operation were transposed (applying the diagonal values in x from the
	// opposite side) for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem );
		bl1_swap_ints( lda, inca );
		bl1_toggle_side( side );
	}

	if ( bl1_is_left( side ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;

			bl1_zewscalv( conj,
			              n_elem,
			              x,       incx,
			              a_begin, inca );
		}
	}
	else
	{
		for ( j = 0; j < n_iter; j++ )
		{
			a_begin = a + j*lda;
			chi     = x + j*incx;
	
			bl1_zscalv( conj,
			            n_elem,
			            chi,
			            a_begin, inca );
		}
	}
}

