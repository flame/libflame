/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sewscalmt( trans1_t trans, integer m, integer n, float* a, integer a_rs, integer a_cs, float* b, integer b_rs, integer b_cs )
{
	float*    a_begin;
	float*    b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_sewscalv( conj,
		              n_elem,
		              a_begin, inca, 
		              b_begin, incb );
	}
}

void bl1_dewscalmt( trans1_t trans, integer m, integer n, double* a, integer a_rs, integer a_cs, double* b, integer b_rs, integer b_cs )
{
	double*   a_begin;
	double*   b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_dewscalv( conj,
		              n_elem,
		              a_begin, inca, 
		              b_begin, incb );
	}
}

void bl1_csewscalmt( trans1_t trans, integer m, integer n, float* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs )
{
	float*    a_begin;
	scomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_csewscalv( conj,
		               n_elem,
		               a_begin, inca, 
		               b_begin, incb );
	}
}

void bl1_cewscalmt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_cewscalv( conj,
		              n_elem,
		              a_begin, inca, 
		              b_begin, incb );
	}
}

void bl1_zdewscalmt( trans1_t trans, integer m, integer n, double* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs )
{
	double*   a_begin;
	dcomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_zdewscalv( conj,
		               n_elem,
		               a_begin, inca, 
		               b_begin, incb );
	}
}

void bl1_zewscalmt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Handle cases where A and B are vectors to ensure that the underlying ewscal
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
		// after a possible transposition, then let's access the matrices by rows
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

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	for ( j = 0; j < n_iter; j++ )
	{
		a_begin = a + j*lda;
		b_begin = b + j*ldb;

		bl1_zewscalv( conj,
		              n_elem,
		              a_begin, inca, 
		              b_begin, incb );
	}
}

