/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sswapmt( trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;
		
		bl1_sswap( n_elem,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bl1_dswapmt( trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;
		
		bl1_dswap( n_elem,
		           a_begin, inca, 
		           b_begin, incb );
	}
}

void bl1_cswapmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;
		
		bl1_cswap( n_elem,
		           a_begin, inca, 
		           b_begin, incb );

		if ( bl1_does_conj( trans ) )
			bl1_cconjv( n_elem,
			            a_begin, inca );

		if ( bl1_does_conj( trans ) )
			bl1_cconjv( n_elem,
			            b_begin, incb );
	}
}

void bl1_zswapmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       j;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying copy
	// gets invoked only once.
	if ( bl1_is_vector( m, n ) )
	{
		// Initialize with values appropriate for vectors.
		n_iter = 1;
		n_elem = bl1_vector_dim( m, n );
		lda    = 1; // multiplied by zero when n_iter == 1; not needed.
		inca   = bl1_vector_inc( trans,             m, n, a_rs, a_cs );
		ldb    = 1; // multiplied by zero when n_iter == 1; not needed.
		incb   = bl1_vector_inc( BLIS1_NO_TRANSPOSE, m, n, b_rs, b_cs );
	}
	else // matrix case
	{
		// Initialize with optimal values for column-major storage.
		n_iter = n;
		n_elem = m;
		lda    = a_cs;
		inca   = a_rs;
		ldb    = b_cs;
		incb   = b_rs;
		
		// Handle the transposition of A.
		if ( bl1_does_trans( trans ) )
		{
			bl1_swap_ints( lda, inca );
		}

		// An optimization: if B is row-major and if A is effectively row-major
		// after a possible transposition, then let's access the matrix by rows
		// instead of by columns for increased spatial locality.
		if ( bl1_is_row_storage( b_rs, b_cs ) )
		{
			if ( ( bl1_is_col_storage( a_rs, a_cs ) && bl1_does_trans( trans ) ) ||
			     ( bl1_is_row_storage( a_rs, a_cs ) && bl1_does_notrans( trans ) ) )
			{
				bl1_swap_ints( n_iter, n_elem );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
			}
		}
	}

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;
		
		bl1_zswap( n_elem,
		           a_begin, inca, 
		           b_begin, incb );

		if ( bl1_does_conj( trans ) )
			bl1_zconjv( n_elem,
			            a_begin, inca );

		if ( bl1_does_conj( trans ) )
			bl1_zconjv( n_elem,
			            b_begin, incb );
	}
}

