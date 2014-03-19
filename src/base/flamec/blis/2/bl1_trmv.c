/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_strmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* a, int a_rs, int a_cs, float* x, int incx )
{
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_screate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
		bl1_toggle_trans( trans );
	}

	bl1_strmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a, lda,
	                x, incx );

	// Free the temporary contiguous matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_dtrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* a, int a_rs, int a_cs, double* x, int incx )
{
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_dcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
		bl1_toggle_trans( trans );
	}

	bl1_dtrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a, lda,
	                x, incx );

	// Free the temporary contiguous matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_ctrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_conj;
	int       incx_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
		bl1_toggle_trans( trans );
	}

	// Initialize with values assuming that trans is not conjnotrans.
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	// Note: strictly speaking, we don't need to create a copy of x since
	// the operation is simpler than, say, gemv. However, we create a copy
	// anyway since in practice it performs better due to increased spatial
	// locality.
	if ( bl1_is_conjnotrans( trans ) )
	{
		x_conj    = bl1_callocv( m );
		incx_conj = 1;

		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bl1_ctrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a,      lda,
	                x_conj, incx_conj );

	// Save the contents of and then free the temporary conjugated x vector.
	if ( bl1_is_conjnotrans( trans ) )
	{
		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    x_conj, incx_conj,
                    x,      incx );

		bl1_cfree( x_conj );
	}

	// Free the temporary contiguous matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_ztrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_conj;
	int       incx_conj;
	int       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( lda, inca );
		bl1_toggle_uplo( uplo );
		bl1_toggle_trans( trans );
	}

	// Initialize with values assuming that trans is not conjnotrans.
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	// Note: strictly speaking, we don't need to create a copy of x since
	// the operation is simpler than, say, gemv. However, we create a copy
	// anyway since in practice it performs better due to increased spatial
	// locality.
	if ( bl1_is_conjnotrans( trans ) )
	{
		x_conj    = bl1_zallocv( m );
		incx_conj = 1;

		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bl1_ztrmv_blas( uplo,
	                trans,
	                diag,
	                m,
	                a,      lda,
	                x_conj, incx_conj );

	// Save the contents of and then free the temporary conjugated x vector.
	if ( bl1_is_conjnotrans( trans ) )
	{
		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    x_conj, incx_conj,
                    x,      incx );

		bl1_zfree( x_conj );
	}

	// Free the temporary contiguous matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_strmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* a, int lda, float* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_strmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_strmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bl1_dtrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* a, int lda, double* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_dtrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_dtrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bl1_ctrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int lda, scomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ctrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ctrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

void bl1_ztrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;
	enum CBLAS_DIAG  cblas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );
	bl1_param_map_to_netlib_diag( diag, &cblas_diag );

	cblas_ztrmv( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             cblas_diag,
	             m,
	             a, lda,
	             x, incx );
#else
	char blas_uplo;
	char blas_trans;
	char blas_diag;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );
	bl1_param_map_to_netlib_diag( diag, &blas_diag );

	F77_ztrmv( &blas_uplo,
	           &blas_trans,
	           &blas_diag,
	           &m,
	           a, &lda,
	           x, &incx );
#endif
}

