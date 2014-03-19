/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_ssyrk( uplo1_t uplo, trans1_t trans, int m, int k, float* alpha, float* a, int a_rs, int a_cs, float* beta, float* c, int c_rs, int c_cs )
{
	uplo1_t    uplo_save = uplo;
	int       m_save    = m;
	float*    a_save    = a;
	float*    c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	int       lda, inca;
	int       ldc, incc;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_screate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_screate_contigmr( uplo,
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
			// requested operation: uplo( C_c ) += A_c * A_c^T
			// effective operation: uplo( C_c ) += A_c * A_c^T
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r^T
			// effective operation: uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c^T
			// effective operation: ~uplo( C_c ) += A_c * A_c^T
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r^T
			// effective operation: ~uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}

	bl1_ssyrk_blas( uplo,
	                trans,
	                m,
	                k,
	                alpha,
	                a, lda,
	                beta,
	                c, ldc );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_sfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_dsyrk( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, double* a, int a_rs, int a_cs, double* beta, double* c, int c_rs, int c_cs )
{
	uplo1_t    uplo_save = uplo;
	int       m_save    = m;
	double*   a_save    = a;
	double*   c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	int       lda, inca;
	int       ldc, incc;

	// Return early if possible.
	if ( bl1_zero_dim2( m, k ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_dcreate_contigmt( trans,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_dcreate_contigmr( uplo,
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
			// requested operation: uplo( C_c ) += A_c * A_c^T
			// effective operation: uplo( C_c ) += A_c * A_c^T
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r^T
			// effective operation: uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c^T
			// effective operation: ~uplo( C_c ) += A_c * A_c^T
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r^T
			// effective operation: ~uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}

	bl1_dsyrk_blas( uplo,
	                trans,
	                m,
	                k,
	                alpha,
	                a, lda,
	                beta,
	                c, ldc );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_dfree_saved_contigmr( uplo_save,
	                          m_save,
	                          m_save,
	                          c_save, c_rs_save, c_cs_save,
	                          &c,     &c_rs,     &c_cs );
}

void bl1_csyrk( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
{
	uplo1_t    uplo_save = uplo;
	int       m_save    = m;
	scomplex* a_save    = a;
	scomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	int       lda, inca;
	int       ldc, incc;

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
			// requested operation: uplo( C_c ) += A_c * A_c^T
			// effective operation: uplo( C_c ) += A_c * A_c^T
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r^T
			// effective operation: uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c^T
			// effective operation: ~uplo( C_c ) += A_c * A_c^T
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r^T
			// effective operation: ~uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}

	bl1_csyrk_blas( uplo,
	                trans,
	                m,
	                k,
	                alpha,
	                a, lda,
	                beta,
	                c, ldc );

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

void bl1_zsyrk( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
{
	uplo1_t    uplo_save = uplo;
	int       m_save    = m;
	dcomplex* a_save    = a;
	dcomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	int       lda, inca;
	int       ldc, incc;
	
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
			// requested operation: uplo( C_c ) += A_c * A_c^T
			// effective operation: uplo( C_c ) += A_c * A_c^T
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation: uplo( C_c ) += A_r * A_r^T
			// effective operation: uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( lda, inca );

			bl1_toggle_trans( trans );
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_c * A_c^T
			// effective operation: ~uplo( C_c ) += A_c * A_c^T
			bl1_swap_ints( ldc, incc );

			bl1_toggle_uplo( uplo );
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			// requested operation:  uplo( C_r ) += A_r * A_r^T
			// effective operation: ~uplo( C_c ) += A_c^T * A_c
			bl1_swap_ints( ldc, incc );
			bl1_swap_ints( lda, inca );

			bl1_toggle_uplo( uplo );
			bl1_toggle_trans( trans );
		}
	}

	bl1_zsyrk_blas( uplo,
	                trans,
	                m,
	                k,
	                alpha,
	                a, lda,
	                beta,
	                c, ldc );

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

void bl1_ssyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, float* alpha, float* a, int lda, float* beta, float* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_ssyrk( cblas_order,
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

	F77_ssyrk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

void bl1_dsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, double* a, int lda, double* beta, double* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_dsyrk( cblas_order,
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

	F77_dsyrk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

void bl1_csyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* beta, scomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_csyrk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             alpha,
	             a, lda,
	             beta,
	             c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_csyrk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

void bl1_zsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* beta, dcomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( trans, &cblas_trans );

	cblas_zsyrk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             alpha,
	             a, lda,
	             beta,
	             c, ldc );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( trans, &blas_trans );

	F77_zsyrk( &blas_uplo,
	           &blas_trans,
	           &m,
	           &k,
	           alpha,
	           a, &lda,
	           beta,
	           c, &ldc );
#endif
}

