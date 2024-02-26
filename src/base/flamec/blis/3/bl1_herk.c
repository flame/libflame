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

void bl1_sherk( uplo1_t uplo, trans1_t trans, integer m, integer k, float* alpha, float* a, integer a_rs, integer a_cs, float* beta, float* c, integer c_rs, integer c_cs )
{
	bl1_ssyrk( uplo,
	           trans,
	           m,
	           k,
	           alpha,
	           a, a_rs, a_cs,
	           beta,
	           c, c_rs, c_cs );
}
void bl1_dherk( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, double* a, integer a_rs, integer a_cs, double* beta, double* c, integer c_rs, integer c_cs )
{
	bl1_dsyrk( uplo,
	           trans,
	           m,
	           k,
	           alpha,
	           a, a_rs, a_cs,
	           beta,
	           c, c_rs, c_cs );
}

void bl1_cherk( uplo1_t uplo, trans1_t trans, integer m, integer k, float* alpha, scomplex* a, integer a_rs, integer a_cs, float* beta, scomplex* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	scomplex* a_save    = a;
	scomplex* c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	float     zero_r = bl1_s0();
	scomplex  one    = bl1_c1();
	scomplex* c_conj;
	integer       lda, inca;
	integer       ldc, incc;
	integer       ldc_conj, incc_conj;
	integer       herk_needs_conj = FALSE;
	
	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_ccreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_c * A_c'
			// effective operation: uplo( C_c ) += A_c * A_c'
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r'
			// effective operation: uplo( C_c ) += conj( A_c' * A_c )
			bl1_swap_ints( lda, inca );

			bl1_toggle_conjtrans( trans );

			herk_needs_conj = TRUE;
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c'
			// effective operation: ~uplo( C_c ) += conj( A_c * A_c' )
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );

			herk_needs_conj = TRUE;
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r'
			// effective operation: ~uplo( C_c ) += A_c' * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_conjtrans( trans );
		}
	}

	// There are two cases where we need to perform the rank-k product and
	// then axpy the result into C with a conjugation. We handle those two
	// cases here.
	if ( herk_needs_conj )
	{
		// We need a temporary matrix for holding the rank-k product.
		c_conj    = bl1_callocm( m, m );
		ldc_conj  = m;
		incc_conj = 1;

		// Compute the rank-k product.
		bl1_cherk_blas( uplo,
		                trans,
		                m,
		                k,
		                alpha,
		                a, lda,
		                &zero_r,
		                c_conj, ldc_conj );

		// Scale C by beta.
		bl1_csscalmr( uplo,
		              m,
		              m,
		              beta,
		              c, incc, ldc );

		// And finally, accumulate the rank-k product in C_conj into C
		// with a conjugation.
		bl1_caxpymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              m,
		              m,
		              &one,
		              c_conj, incc_conj, ldc_conj,
		              c,      incc,      ldc );

		// Free the temporary matrix for C.
		bl1_cfree( c_conj );
	}
	else
	{
		bl1_cherk_blas( uplo,
		                trans,
		                m,
		                k,
		                alpha,
		                a, lda,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_cfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_zherk( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, dcomplex* a, integer a_rs, integer a_cs, double* beta, dcomplex* c, integer c_rs, integer c_cs )
{
	uplo1_t    uplo_save = uplo;
	integer       m_save    = m;
	dcomplex* a_save    = a;
	dcomplex* c_save    = c;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       c_rs_save = c_rs;
	integer       c_cs_save = c_cs;
	double    zero_r = bl1_d0();
	dcomplex  one    = bl1_z1();
	dcomplex* c_conj;
	integer       lda, inca;
	integer       ldc, incc;
	integer       ldc_conj, incc_conj;
	integer       herk_needs_conj = FALSE;
	
	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_zcreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      c_save, c_rs_save, c_cs_save,
	                      &c,     &c_rs,     &c_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldc  = c_cs;
	incc = c_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_c * A_c'
			// effective operation: uplo( C_c ) += A_c * A_c'
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r'
			// effective operation: uplo( C_c ) += conj( A_c' * A_c )
			bl1_swap_ints( lda, inca );

			bl1_toggle_conjtrans( trans );

			herk_needs_conj = TRUE;
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c'
			// effective operation: ~uplo( C_c ) += conj( A_c * A_c' )
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );

			herk_needs_conj = TRUE;
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r'
			// effective operation: ~uplo( C_c ) += A_c' * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_conjtrans( trans );
		}
	}

	// There are two cases where we need to perform the rank-k product and
	// then axpy the result into C with a conjugation. We handle those two
	// cases here.
	if ( herk_needs_conj )
	{
		// We need a temporary matrix for holding the rank-k product.
		c_conj    = bl1_zallocm( m, m );
		ldc_conj  = m;
		incc_conj = 1;

		// Compute the rank-k product.
		bl1_zherk_blas( uplo,
		                trans,
		                m,
		                k,
		                alpha,
		                a, lda,
		                &zero_r,
		                c_conj, ldc_conj );

		// Scale C by beta.
		bl1_zdscalmr( uplo,
		              m,
		              m,
		              beta,
		              c, incc, ldc );
		
		// And finally, accumulate the rank-k product in C_conj into C
		// with a conjugation.
		bl1_zaxpymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              m,
		              m,
		              &one,
		              c_conj, incc_conj, ldc_conj,
		              c,      incc,      ldc );

		// Free the temporary matrix for C.
		bl1_zfree( c_conj );
	}
	else
	{
		bl1_zherk_blas( uplo,
		                trans,
		                m,
		                k,
		                alpha,
		                a, lda,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_zfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

// --- Classic routine wrappers ---

void bl1_cherk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, float* alpha, scomplex* a, integer lda, float* beta, scomplex* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_cherk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             *alpha,
	             a, lda,
	             *beta,
	             c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_cherk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

void bl1_zherk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, dcomplex* a, integer lda, double* beta, dcomplex* c, integer ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_zherk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             *alpha,
	             a, lda,
	             *beta,
	             c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_zherk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

