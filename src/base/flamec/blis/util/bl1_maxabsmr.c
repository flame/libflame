/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_smaxabsmr( uplo1_t uplo, integer m, integer n, float* a, integer a_rs, integer a_cs, float* maxabs )
{
	float     zero = bl1_s0();
	float*    a_begin;
	float     maxabs_cand;
	float     maxabs_temp;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) { *maxabs = zero; return; }

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize the maximum absolute value candidate to the first element.
	bl1_sabsval2( a, &maxabs_cand );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_smaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			bl1_smaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}

	*maxabs = maxabs_cand;
}

void bl1_dmaxabsmr( uplo1_t uplo, integer m, integer n, double* a, integer a_rs, integer a_cs, double* maxabs )
{
	double    zero = bl1_d0();
	double*   a_begin;
	double    maxabs_cand;
	double    maxabs_temp;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) { *maxabs = zero; return; }

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize the maximum absolute value candidate to the first element.
	bl1_dabsval2( a, &maxabs_cand );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_dmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			bl1_dmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}

	*maxabs = maxabs_cand;
}

void bl1_cmaxabsmr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float* maxabs )
{
	float     zero = bl1_d0();
	scomplex* a_begin;
	float     maxabs_cand;
	float     maxabs_temp;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) { *maxabs = zero; return; }

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize the maximum absolute value candidate to the first element.
	bl1_csabsval2( a, &maxabs_cand );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_cmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			bl1_cmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}

	*maxabs = maxabs_cand;
}

void bl1_zmaxabsmr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double* maxabs )
{
	double    zero = bl1_d0();
	dcomplex* a_begin;
	double    maxabs_cand;
	double    maxabs_temp;
	integer       inca, lda;
	integer       n_iter;
	integer       n_elem_max;
	integer       n_elem;
	integer       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) { *maxabs = zero; return; }

	// Initialize with optimal values for column-major storage.
	n_iter     = n;
	n_elem_max = m;
	lda        = a_cs;
	inca       = a_rs;

	// An optimization: if A is row-major, then let's access the matrix by
	// rows instead of by columns for increased spatial locality.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( n_iter, n_elem_max );
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
	}

	// Initialize the maximum absolute value candidate to the first element.
	bl1_zdabsval2( a, &maxabs_cand );

	if ( bl1_is_upper( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;

			bl1_zmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}
	else // if ( bl1_is_lower( uplo ) )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_max( 0, n_elem_max - j );
			a_begin = a + j*lda + j*inca;

			bl1_zmaxabsv( n_elem,
			              a_begin, inca,
			              &maxabs_temp );

			if ( maxabs_temp > maxabs_cand ) maxabs_cand = maxabs_temp;
		}
	}

	*maxabs = maxabs_cand;
}

