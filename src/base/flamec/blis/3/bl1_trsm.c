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

void bl1_strsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float* alpha, float* a, integer a_rs, integer a_cs, float* b, integer b_rs, integer b_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	float*    a_save    = a;
	float*    b_save    = b;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       dim_a;
	integer       lda, inca;
	integer       ldb, incb;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_screate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_screate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr( uplo( A_c ) ) * B_c
			// effective operation: B_c := tr( uplo( A_c ) ) * B_c
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr(  uplo( A_r ) )   * B_c
			// effective operation: B_c := tr( ~uplo( A_c ) )^T * B_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_c ) ) * B_r
			// effective operation: B_c := B_c * tr( uplo( A_c ) )^T
			bl1_swap_ints( ldb, incb );

			bl1_swap_ints( m, n );

			bl1_toggle_side( side );
			bl1_toggle_trans( trans );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_r ) ) * B_r
			// effective operation: B_c := B_c * tr( ~uplo( A_c ) )
			bl1_swap_ints( ldb, incb );
			bl1_swap_ints( lda, inca );

			bl1_swap_ints( m, n );

			bl1_toggle_uplo( uplo );
			bl1_toggle_side( side );
		}
	}

	bl1_strsm_blas( side,
	                uplo,
	                trans,
	                diag,
	                m,
	                n,
	                alpha,
	                a, lda,
	                b, ldb );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_sfree_saved_contigm( m_save,
	                         n_save,
	                         b_save, b_rs_save, b_cs_save,
	                         &b,     &b_rs,     &b_cs );
}

void bl1_dtrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double* alpha, double* a, integer a_rs, integer a_cs, double* b, integer b_rs, integer b_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	double*   a_save    = a;
	double*   b_save    = b;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	integer       dim_a;
	integer       lda, inca;
	integer       ldb, incb;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_dcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_dcreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr( uplo( A_c ) ) * B_c
			// effective operation: B_c := tr( uplo( A_c ) ) * B_c
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr(  uplo( A_r ) )   * B_c
			// effective operation: B_c := tr( ~uplo( A_c ) )^T * B_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_c ) ) * B_r
			// effective operation: B_c := B_c * tr( uplo( A_c ) )^T
			bl1_swap_ints( ldb, incb );

			bl1_swap_ints( m, n );

			bl1_toggle_side( side );
			bl1_toggle_trans( trans );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_r ) ) * B_r
			// effective operation: B_c := B_c * tr( ~uplo( A_c ) )
			bl1_swap_ints( ldb, incb );
			bl1_swap_ints( lda, inca );

			bl1_swap_ints( m, n );

			bl1_toggle_uplo( uplo );
			bl1_toggle_side( side );
		}
	}

	bl1_dtrsm_blas( side,
	                uplo,
	                trans,
	                diag,
	                m,
	                n,
	                alpha,
	                a, lda,
	                b, ldb );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_dfree_saved_contigm( m_save,
	                         n_save,
	                         b_save, b_rs_save, b_cs_save,
	                         &b,     &b_rs,     &b_cs );
}

void bl1_ctrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	scomplex* a_save    = a;
	scomplex* b_save    = b;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	scomplex* a_conj;
	integer       dim_a;
	integer       lda, inca;
	integer       ldb, incb;
	integer       lda_conj, inca_conj;
	integer       a_was_copied;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_ccreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_ccreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	// Figure out whether A was copied to contiguous memory. This is used to
	// prevent redundant copying.
	a_was_copied = ( a != a_save );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr( uplo( A_c ) ) * B_c
			// effective operation: B_c := tr( uplo( A_c ) ) * B_c
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr(  uplo( A_r ) )   * B_c
			// effective operation: B_c := tr( ~uplo( A_c ) )^T * B_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_c ) ) * B_r
			// effective operation: B_c := B_c * tr( uplo( A_c ) )^T
			bl1_swap_ints( ldb, incb );

			bl1_swap_ints( m, n );

			bl1_toggle_side( side );
			bl1_toggle_trans( trans );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_r ) ) * B_r
			// effective operation: B_c := B_c * tr( ~uplo( A_c ) )
			bl1_swap_ints( ldb, incb );
			bl1_swap_ints( lda, inca );

			bl1_swap_ints( m, n );

			bl1_toggle_uplo( uplo );
			bl1_toggle_side( side );
		}
	}

	// Initialize with values assuming that trans is not conjnotrans.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// We want to handle the conjnotrans case. The easiest way to do so is
	// by making a conjugated copy of A.
	if ( bl1_is_conjnotrans( trans ) && !a_was_copied )
	{
		integer dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		a_conj    = bl1_callocm( dim_a, dim_a );
		lda_conj  = dim_a;
		inca_conj = 1;

		bl1_ccopymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              dim_a,
		              dim_a,
		              a,      inca,      lda,
		              a_conj, inca_conj, lda_conj );
	}
	else if ( bl1_is_conjnotrans( trans ) && a_was_copied )
	{
		integer dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		bl1_cconjmr( uplo,
		             dim_a,
		             dim_a,
		             a_conj, inca_conj, lda_conj );
	}


	bl1_ctrsm_blas( side,
	                uplo,
	                trans,
	                diag,
	                m,
	                n,
	                alpha,
	                a_conj, lda_conj,
	                b,      ldb );

	if ( bl1_is_conjnotrans( trans ) && !a_was_copied )
		bl1_cfree( a_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_cfree_saved_contigm( m_save,
	                         n_save,
	                         b_save, b_rs_save, b_cs_save,
	                         &b,     &b_rs,     &b_cs );
}

void bl1_ztrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	dcomplex* a_save    = a;
	dcomplex* b_save    = b;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       b_rs_save = b_rs;
	integer       b_cs_save = b_cs;
	dcomplex* a_conj;
	integer       dim_a;
	integer       lda, inca;
	integer       ldb, incb;
	integer       lda_conj, inca_conj;
	integer       a_was_copied;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_zcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_zcreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	// Figure out whether A was copied to contiguous memory. This is used to
	// prevent redundant copying.
	a_was_copied = ( a != a_save );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;

	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr( uplo( A_c ) ) * B_c
			// effective operation: B_c := tr( uplo( A_c ) ) * B_c
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_c := tr(  uplo( A_r ) )   * B_c
			// effective operation: B_c := tr( ~uplo( A_c ) )^T * B_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_c ) ) * B_r
			// effective operation: B_c := B_c * tr( uplo( A_c ) )^T
			bl1_swap_ints( ldb, incb );

			bl1_swap_ints( m, n );

			bl1_toggle_side( side );
			bl1_toggle_trans( trans );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: B_r := tr( uplo( A_r ) ) * B_r
			// effective operation: B_c := B_c * tr( ~uplo( A_c ) )
			bl1_swap_ints( ldb, incb );
			bl1_swap_ints( lda, inca );

			bl1_swap_ints( m, n );

			bl1_toggle_side( side );
			bl1_toggle_uplo( uplo );
		}
	}

	// Initialize with values assuming that trans is not conjnotrans.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// We want to handle the conjnotrans case. The easiest way to do so is
	// by making a conjugated copy of A.
	if ( bl1_is_conjnotrans( trans ) && !a_was_copied )
	{
		integer dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		a_conj    = bl1_zallocm( dim_a, dim_a );
		lda_conj  = dim_a;
		inca_conj = 1;

		bl1_zcopymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              dim_a,
		              dim_a,
		              a,      inca,      lda,
		              a_conj, inca_conj, lda_conj );
	}
	else if ( bl1_is_conjnotrans( trans ) && a_was_copied )
	{
		integer dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		bl1_zconjmr( uplo,
		             dim_a,
		             dim_a,
		             a_conj, inca_conj, lda_conj );
	}

	bl1_ztrsm_blas( side,
	                uplo,
	                trans,
	                diag,
	                m,
	                n,
	                alpha,
	                a_conj, lda_conj,
	                b,      ldb );

	if ( bl1_is_conjnotrans( trans ) && !a_was_copied )
		bl1_zfree( a_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_zfree_saved_contigm( m_save,
	                         n_save,
	                         b_save, b_rs_save, b_cs_save,
	                         &b,     &b_rs,     &b_cs );
}

// --- Classic routine wrappers ---

void bl1_strsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float* alpha, float* a, integer lda, float* b, integer ldb )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_strsm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             b, ldb );
#else
	char blas_side;
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_strsm( &blas_side,
	           &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb );
#endif
}

void bl1_dtrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double* alpha, double* a, integer lda, double* b, integer ldb )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_dtrsm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             b, ldb );
#else
	char blas_side;
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_dtrsm( &blas_side,
	           &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb );
#endif
}

void bl1_ctrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ctrsm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             n,
	             alpha,
	             a, lda,
	             b, ldb );
#else
	char blas_side;
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ctrsm( &blas_side,
	           &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb );
#endif
}

void bl1_ztrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ztrsm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             n,
	             alpha,
	             a, lda,
	             b, ldb );
#else
	char blas_side;
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ztrsm( &blas_side,
	           &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb );
#endif
}

