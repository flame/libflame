/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_ssetmr( uplo1_t uplo, integer m, integer n, float* sigma, float* a, integer a_rs, integer a_cs )
{
	float*    a_begin;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			bl1_ssetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bl1_ssetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bl1_dsetmr( uplo1_t uplo, integer m, integer n, double* sigma, double* a, integer a_rs, integer a_cs )
{
	double*   a_begin;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			bl1_dsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bl1_dsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bl1_csetmr( uplo1_t uplo, integer m, integer n, scomplex* sigma, scomplex* a, integer a_rs, integer a_cs )
{
	scomplex* a_begin;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			bl1_csetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bl1_csetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

void bl1_zsetmr( uplo1_t uplo, integer m, integer n, dcomplex* sigma, dcomplex* a, integer a_rs, integer a_cs )
{
	dcomplex* a_begin;
	integer       lda, inca;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns to increase spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}
	
	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j, n_elem_max );
			a_begin = a + j*lda;

			bl1_zsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j - 1 );
			a_begin = a + j*lda + (j + 1)*inca;

			bl1_zsetv( n_elem,
			           sigma,
			           a_begin, inca );
		}
	}
}

