/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_ssyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, float* alpha, float* a, integer a_rs, integer a_cs, float* b, integer b_rs, integer b_cs, float* beta, float* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	float*    a_save    = a;
	float*    b_save    = b;
	float*    c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	float*    a_copy;
	float*    b_copy;
	integer       lda, inca;
	integer       ldb, incb;
	integer       ldc, incc;
	integer       lda_copy, inca_copy;
	integer       ldb_copy, incb_copy;
	integer       syr2k_needs_copya = FALSE;
	integer       syr2k_needs_copyb = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_screate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_screate_contigmt( trans,
	                      m,
	                      k,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_screate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_r' + B_r * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_c' + B_c * A_r'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copya = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_r' + B_r * A_r'
				// requested operation: uplo( C_c ) += conj( A_c' * B_c + B_c' * A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( trans );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_c' + B_c * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_r' + B_r * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copyb = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_c' + B_c * A_r'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copya = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_r' + B_r * A_r'
				// requested operation: ~uplo( C_c ) += A_c' * B_c + B_c' * A_c
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_uplo( uplo );
				bl1_toggle_trans( trans );
			}
		}
	}

	a_copy    = a;
	lda_copy  = lda;
	inca_copy = inca;
	
	// There are two cases where we need to copy A column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copya )
	{
		integer m_a;
		integer n_a;

		// Determine the dimensions of A according to the value of trans. We
		// need this in order to set the leading dimension of the copy of A.
		bl1_set_dims_with_trans( trans, m, k, &m_a, &n_a );

		// We need a temporary matrix to hold a column-major copy of A.
		a_copy    = bl1_sallocm( m, k );
		lda_copy  = m_a;
		inca_copy = 1;

		// Copy the contents of A into A_copy.
		bl1_scopymt( BLIS1_NO_TRANSPOSE,
                     m_a,
                     n_a,
		             a,      inca,      lda,
		             a_copy, inca_copy, lda_copy );
	}
	
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;

	// There are two cases where we need to copy B column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copyb )
	{
		integer m_b;
		integer n_b;

		// Determine the dimensions of B according to the value of trans. We
		// need this in order to set the leading dimension of the copy of B.
		bl1_set_dims_with_trans( trans, m, k, &m_b, &n_b );

		// We need a temporary matrix to hold a column-major copy of B.
		b_copy    = bl1_sallocm( m, k );
		ldb_copy  = m_b;
		incb_copy = 1;

		// Copy the contents of B into B_copy.
		bl1_scopymt( BLIS1_NO_TRANSPOSE,
                     m_b,
                     n_b,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	bl1_ssyr2k_blas( uplo,
	                 trans,
	                 m,
	                 k,
	                 alpha,
	                 a_copy, lda_copy,
	                 b_copy, ldb_copy,
	                 beta,
	                 c, ldc );

	if ( syr2k_needs_copya )
		bl1_sfree( a_copy );

	if ( syr2k_needs_copyb )
		bl1_sfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_sfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_sfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_dsyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, double* a, integer a_rs, integer a_cs, double* b, integer b_rs, integer b_cs, double* beta, double* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	double*   a_save    = a;
	double*   b_save    = b;
	double*   c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	double*   a_copy;
	double*   b_copy;
	integer       lda, inca;
	integer       ldb, incb;
	integer       ldc, incc;
	integer       lda_copy, inca_copy;
	integer       ldb_copy, incb_copy;
	integer       syr2k_needs_copya = FALSE;
	integer       syr2k_needs_copyb = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_dcreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_dcreate_contigmt( trans,
	                      m,
	                      k,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_dcreate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_r' + B_r * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_c' + B_c * A_r'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copya = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_r' + B_r * A_r'
				// requested operation: uplo( C_c ) += conj( A_c' * B_c + B_c' * A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( trans );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_c' + B_c * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_r' + B_r * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copyb = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_c' + B_c * A_r'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copya = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_r' + B_r * A_r'
				// requested operation: ~uplo( C_c ) += A_c' * B_c + B_c' * A_c
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_uplo( uplo );
				bl1_toggle_trans( trans );
			}
		}
	}

	a_copy    = a;
	lda_copy  = lda;
	inca_copy = inca;
	
	// There are two cases where we need to copy A column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copya )
	{
		integer m_a;
		integer n_a;

		// Determine the dimensions of A according to the value of trans. We
		// need this in order to set the leading dimension of the copy of A.
		bl1_set_dims_with_trans( trans, m, k, &m_a, &n_a );

		// We need a temporary matrix to hold a column-major copy of A.
		a_copy    = bl1_dallocm( m, k );
		lda_copy  = m_a;
		inca_copy = 1;

		// Copy the contents of A into A_copy.
		bl1_dcopymt( BLIS1_NO_TRANSPOSE,
                     m_a,
                     n_a,
		             a,      inca,      lda,
		             a_copy, inca_copy, lda_copy );
	}
	
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;

	// There are two cases where we need to copy B column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copyb )
	{
		integer m_b;
		integer n_b;

		// Determine the dimensions of B according to the value of trans. We
		// need this in order to set the leading dimension of the copy of B.
		bl1_set_dims_with_trans( trans, m, k, &m_b, &n_b );

		// We need a temporary matrix to hold a column-major copy of B.
		b_copy    = bl1_dallocm( m, k );
		ldb_copy  = m_b;
		incb_copy = 1;

		// Copy the contents of B into B_copy.
		bl1_dcopymt( BLIS1_NO_TRANSPOSE,
                     m_b,
                     n_b,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	bl1_dsyr2k_blas( uplo,
	                 trans,
	                 m,
	                 k,
	                 alpha,
	                 a_copy, lda_copy,
	                 b_copy, ldb_copy,
	                 beta,
	                 c, ldc );

	if ( syr2k_needs_copya )
		bl1_dfree( a_copy );

	if ( syr2k_needs_copyb )
		bl1_dfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_dfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_dfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_csyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	scomplex* a_save    = a;
	scomplex* b_save    = b;
	scomplex* c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	scomplex* a_copy;
	scomplex* b_copy;
	integer       lda, inca;
	integer       ldb, incb;
	integer       ldc, incc;
	integer       lda_copy, inca_copy;
	integer       ldb_copy, incb_copy;
	integer       syr2k_needs_copya = FALSE;
	integer       syr2k_needs_copyb = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_ccreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_ccreate_contigmt( trans,
	                      m,
	                      k,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_r' + B_r * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_c' + B_c * A_r'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copya = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_r' + B_r * A_r'
				// requested operation: uplo( C_c ) += conj( A_c' * B_c + B_c' * A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( trans );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_c' + B_c * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_r' + B_r * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copyb = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_c' + B_c * A_r'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copya = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_r' + B_r * A_r'
				// requested operation: ~uplo( C_c ) += A_c' * B_c + B_c' * A_c
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_uplo( uplo );
				bl1_toggle_trans( trans );
			}
		}
	}

	a_copy    = a;
	lda_copy  = lda;
	inca_copy = inca;
	
	// There are two cases where we need to copy A column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copya )
	{
		integer m_a;
		integer n_a;

		// Determine the dimensions of A according to the value of trans. We
		// need this in order to set the leading dimension of the copy of A.
		bl1_set_dims_with_trans( trans, m, k, &m_a, &n_a );

		// We need a temporary matrix to hold a column-major copy of A.
		a_copy    = bl1_callocm( m, k );
		lda_copy  = m_a;
		inca_copy = 1;

		// Copy the contents of A into A_copy.
		bl1_ccopymt( BLIS1_NO_TRANSPOSE,
                     m_a,
                     n_a,
		             a,      inca,      lda,
		             a_copy, inca_copy, lda_copy );
	}
	
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;

	// There are two cases where we need to copy B column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copyb )
	{
		integer m_b;
		integer n_b;

		// Determine the dimensions of B according to the value of trans. We
		// need this in order to set the leading dimension of the copy of B.
		bl1_set_dims_with_trans( trans, m, k, &m_b, &n_b );

		// We need a temporary matrix to hold a column-major copy of B.
		b_copy    = bl1_callocm( m, k );
		ldb_copy  = m_b;
		incb_copy = 1;

		// Copy the contents of B into B_copy.
		bl1_ccopymt( BLIS1_NO_TRANSPOSE,
                     m_b,
                     n_b,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	bl1_csyr2k_blas( uplo,
	                 trans,
	                 m,
	                 k,
	                 alpha,
	                 a_copy, lda_copy,
	                 b_copy, ldb_copy,
	                 beta,
	                 c, ldc );

	if ( syr2k_needs_copya )
		bl1_cfree( a_copy );

	if ( syr2k_needs_copyb )
		bl1_cfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_cfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_cfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_zsyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	dcomplex* a_save    = a;
	dcomplex* b_save    = b;
	dcomplex* c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	dcomplex* a_copy;
	dcomplex* b_copy;
	integer       lda, inca;
	integer       ldb, incb;
	integer       ldc, incc;
	integer       lda_copy, inca_copy;
	integer       ldb_copy, incb_copy;
	integer       syr2k_needs_copya = FALSE;
	integer       syr2k_needs_copyb = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_zcreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_zcreate_contigmt( trans,
	                      m,
	                      k,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_c * B_r' + B_r * A_c'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_c' + B_c * A_r'
				// requested operation: uplo( C_c ) += A_c * B_c' + B_c * A_c'
				syr2k_needs_copya = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: uplo( C_c ) += A_r * B_r' + B_r * A_r'
				// requested operation: uplo( C_c ) += conj( A_c' * B_c + B_c' * A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( trans );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_c' + B_c * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_c * B_r' + B_r * A_c'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copyb = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_c' + B_c * A_r'
				// requested operation: ~uplo( C_c ) += conj( A_c * B_c' + B_c * A_c' )
				syr2k_needs_copya = TRUE;

				bl1_swap_ints( ldc, incc );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation:  uplo( C_r ) += A_r * B_r' + B_r * A_r'
				// requested operation: ~uplo( C_c ) += A_c' * B_c + B_c' * A_c
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_uplo( uplo );
				bl1_toggle_trans( trans );
			}
		}
	}

	a_copy    = a;
	lda_copy  = lda;
	inca_copy = inca;
	
	// There are two cases where we need to copy A column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copya )
	{
		integer m_a;
		integer n_a;

		// Determine the dimensions of A according to the value of trans. We
		// need this in order to set the leading dimension of the copy of A.
		bl1_set_dims_with_trans( trans, m, k, &m_a, &n_a );

		// We need a temporary matrix to hold a column-major copy of A.
		a_copy    = bl1_zallocm( m, k );
		lda_copy  = m_a;
		inca_copy = 1;

		// Copy the contents of A into A_copy.
		bl1_zcopymt( BLIS1_NO_TRANSPOSE,
                     m_a,
                     n_a,
		             a,      inca,      lda,
		             a_copy, inca_copy, lda_copy );
	}
	
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;

	// There are two cases where we need to copy B column-major storage.
	// We handle those two cases here.
	if ( syr2k_needs_copyb )
	{
		integer m_b;
		integer n_b;

		// Determine the dimensions of B according to the value of trans. We
		// need this in order to set the leading dimension of the copy of B.
		bl1_set_dims_with_trans( trans, m, k, &m_b, &n_b );

		// We need a temporary matrix to hold a column-major copy of B.
		b_copy    = bl1_zallocm( m, k );
		ldb_copy  = m_b;
		incb_copy = 1;

		// Copy the contents of B into B_copy.
		bl1_zcopymt( BLIS1_NO_TRANSPOSE,
                     m_b,
                     n_b,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	bl1_zsyr2k_blas( uplo,
	                 trans,
	                 m,
	                 k,
	                 alpha,
	                 a_copy, lda_copy,
	                 b_copy, ldb_copy,
	                 beta,
	                 c, ldc );

	if ( syr2k_needs_copya )
		bl1_zfree( a_copy );

	if ( syr2k_needs_copyb )
		bl1_zfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_zfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_zfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

// --- Classic routine wrappers ---

void bl1_ssyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, float* alpha, float* a, integer lda, float* b, integer ldb, float* beta, float* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_ssyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              *alpha,
	              a, lda,
	              b, ldb,
	              *beta,
	              c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_ssyr2k( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            a, &lda,
	            b, &ldb,
	            beta,
	            c, &ldc );
#endif
}

void bl1_dsyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, double* a, integer lda, double* b, integer ldb, double* beta, double* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_dsyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              *alpha,
	              a, lda,
	              b, ldb,
	              *beta,
	              c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_dsyr2k( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            a, &lda,
	            b, &ldb,
	            beta,
	            c, &ldc );
#endif
}

void bl1_csyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, scomplex* beta, scomplex* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_csyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              a, lda,
	              b, ldb,
	              beta,
	              c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_csyr2k( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            a, &lda,
	            b, &ldb,
	            beta,
	            c, &ldc );
#endif
}

void bl1_zsyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, dcomplex* beta, dcomplex* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_zsyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              a, lda,
	              b, ldb,
	              beta,
	              c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	// BLAS doesn't recognize the conjugate-transposition constant for syr2k,
	// so we have to map it down to regular transposition.
	if ( bl1_is_conjtrans( trans ) ) trans = BLIS1_TRANSPOSE;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_zsyr2k( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            a, &lda,
	            b, &ldb,
	            beta,
	            c, &ldc );
#endif
}

