/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sher( uplo1_t uplo, conj1_t conj, int m, float* alpha, float* x, int incx, float* a, int a_rs, int a_cs )
{
	bl1_ssyr( uplo,
	          m,
	          alpha,
	          x, incx,
	          a, a_rs, a_cs );
}

void bl1_dher( uplo1_t uplo, conj1_t conj, int m, double* alpha, double* x, int incx, double* a, int a_rs, int a_cs )
{
	bl1_dsyr( uplo,
	          m,
	          alpha,
	          x, incx,
	          a, a_rs, a_cs );
}

void bl1_cher( uplo1_t uplo, conj1_t conj, int m, float* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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
		bl1_toggle_conj( conj );
	}

	// Initialize with values assuming no conjugation of ( x * x' ).
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the case where ( x * x' ) is conjugated, but
	// without explicitly conjugating the matrix. To do so, we leverage
	// the fact that computing the product conj( x * x' ) is equivalent
	// to computing ( conj(x) * conj(x)' ), since ( x * x' ) is Hermitian.
	if ( bl1_is_conj( conj ) )
	{
		x_conj    = bl1_callocv( m );
		incx_conj = 1;

		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bl1_cher_blas( uplo,
	               m,
	               alpha,
	               x_conj, incx_conj,
	               a,      lda );

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conj ) )
		bl1_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_zher( uplo1_t uplo, conj1_t conj, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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
		bl1_toggle_conj( conj );
	}

	// Initialize with values assuming no conjugation of ( x * x' ).
	x_conj    = x;
	incx_conj = incx;

	// We want to handle the case where ( x * x' ) is conjugated, but
	// without explicitly conjugating the matrix. To do so, we leverage
	// the fact that computing the product conj( x * x' ) is equivalent
	// to computing ( conj(x) * conj(x)' ), since ( x * x' ) is Hermitian.
	if ( bl1_is_conj( conj ) )
	{
		x_conj    = bl1_zallocv( m );
		incx_conj = 1;

		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	bl1_zher_blas( uplo,
	               m,
	               alpha,
	               x_conj, incx_conj,
	               a,      lda );

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conj ) )
		bl1_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_cher_blas( uplo1_t uplo, int m, float* alpha, scomplex* x, int incx, scomplex* a, int lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_cher( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_cher( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

void bl1_zher_blas( uplo1_t uplo, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zher( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zher( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

