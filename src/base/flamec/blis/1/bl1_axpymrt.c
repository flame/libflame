/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_saxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float* alpha, float* a, integer a_rs, integer a_cs, float* b, integer b_rs, integer b_cs )
{
	float*    a_begin;
	float*    b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       n_elem_max;
	integer       n_elem_is_descending;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bl1_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bl1_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bl1_does_trans( trans ) )
	{
		bl1_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bl1_saxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bl1_saxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_daxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double* alpha, double* a, integer a_rs, integer a_cs, double* b, integer b_rs, integer b_cs )
{
	double*   a_begin;
	double*   b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       n_elem_max;
	integer       n_elem_is_descending;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bl1_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bl1_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bl1_does_trans( trans ) )
	{
		bl1_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bl1_daxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bl1_daxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_caxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       n_elem_max;
	integer       n_elem_is_descending;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bl1_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bl1_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bl1_does_trans( trans ) )
	{
		bl1_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bl1_caxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bl1_caxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_zaxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	integer       lda, inca;
	integer       ldb, incb;
	integer       n_iter;
	integer       n_elem;
	integer       n_elem_max;
	integer       n_elem_is_descending;
	integer       j;
	conj1_t    conj;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// Initialize variables based on storage format of B and value of uplo.
	if      ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = m;
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = TRUE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = n;
			n_elem_max = bl1_min( m, n );
			lda        = a_cs;
			inca       = a_rs;
			ldb        = b_cs;
			incb       = b_rs;
			n_elem_is_descending = FALSE;
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_lower( uplo ) )
		{
			n_iter     = m;
			n_elem_max = bl1_min( m, n );
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = FALSE;
		}
		else // if ( bl1_is_upper( uplo ) )
		{
			n_iter     = bl1_min( m, n );
			n_elem_max = n;
			lda        = a_rs;
			inca       = a_cs;
			ldb        = b_rs;
			incb       = b_cs;
			n_elem_is_descending = TRUE;
		}
	}

	// Swap lda and inca if we're doing a transpose.
	if ( bl1_does_trans( trans ) )
	{
		bl1_swap_ints( lda, inca );
	}

	// Extract conj component from trans parameter.
	conj = bl1_proj_trans1_to_conj( trans );

	// Choose the loop based on whether n_elem will be shrinking or growing
	// with each iteration.
	if ( n_elem_is_descending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = n_elem_max - j;
			a_begin = a + j*lda + j*inca;
			b_begin = b + j*ldb + j*incb;
		
			bl1_zaxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
	else // if ( n_elem_is_ascending )
	{
		for ( j = 0; j < n_iter; j++ )
		{
			n_elem  = bl1_min( j + 1, n_elem_max );
			a_begin = a + j*lda;
			b_begin = b + j*ldb;
		
			bl1_zaxpyv( conj,
			            n_elem,
			            alpha,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

