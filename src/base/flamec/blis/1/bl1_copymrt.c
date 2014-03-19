/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_scopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_scopyv( conj,
			            n_elem,
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
		
			bl1_scopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_dcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_dcopyv( conj,
			            n_elem,
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
		
			bl1_dcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_ccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_ccopyv( conj,
			            n_elem,
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
		
			bl1_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

void bl1_zcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_zcopyv( conj,
			            n_elem,
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
		
			bl1_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// --- Mixed-datatype and general stride copy routines---------------

// ss
void bl1_sscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	float*    a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_scopyv( conj,
			            n_elem,
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
		
			bl1_scopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// sd
void bl1_sdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	float*    a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_sdcopyv( conj,
			             n_elem,
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
		
			bl1_sdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// sc
void bl1_sccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_sccopyv( conj,
			             n_elem,
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
		
			bl1_sccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// sz
void bl1_szcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	float*    a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_szcopyv( conj,
			             n_elem,
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
		
			bl1_szcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// ds
void bl1_dscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	double*   a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_dscopyv( conj,
			             n_elem,
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
		
			bl1_dscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// dd
void bl1_ddcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	double*   a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_dcopyv( conj,
			            n_elem,
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
		
			bl1_dcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// dc
void bl1_dccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_dccopyv( conj,
			             n_elem,
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
		
			bl1_dccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// dz
void bl1_dzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	double*   a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_dzcopyv( conj,
			             n_elem,
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
		
			bl1_dzcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cs
void bl1_cscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_cscopyv( conj,
			             n_elem,
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
		
			bl1_cscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cd
void bl1_cdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_cdcopyv( conj,
			             n_elem,
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
		
			bl1_cdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// cc
void bl1_cccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_ccopyv( conj,
			            n_elem,
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
		
			bl1_ccopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

// cz
void bl1_czcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	scomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_czcopyv( conj,
			             n_elem,
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
		
			bl1_czcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zs
void bl1_zscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	float*    b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_zscopyv( conj,
			             n_elem,
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
		
			bl1_zscopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zd
void bl1_zdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	double*   b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_zdcopyv( conj,
			             n_elem,
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
		
			bl1_zdcopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zc
void bl1_zccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	scomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_zccopyv( conj,
			             n_elem,
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
		
			bl1_zccopyv( conj,
			             n_elem,
			             a_begin, inca,
			             b_begin, incb );
		}
	}
}

// zz
void bl1_zzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs )
{
	dcomplex* a_begin;
	dcomplex* b_begin;
	int       lda, inca;
	int       ldb, incb;
	int       n_iter;
	int       n_elem;
	int       n_elem_max;
	int       n_elem_is_descending;
	int       j;
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
		
			bl1_zcopyv( conj,
			            n_elem,
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
		
			bl1_zcopyv( conj,
			            n_elem,
			            a_begin, inca,
			            b_begin, incb );
		}
	}
}

