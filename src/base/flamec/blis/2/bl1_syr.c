
#include "blis1.h"

void bl1_ssyr( uplo1_t uplo, int m, float* alpha, float* x, int incx, float* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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
	}

	bl1_ssyr_blas( uplo,
	               m,
	               alpha,
	               x, incx,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_sfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_dsyr( uplo1_t uplo, int m, double* alpha, double* x, int incx, double* a, int a_rs, int a_cs )
{
	int       m_save    = m;
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
	}

	bl1_dsyr_blas( uplo,
	               m,
	               alpha,
	               x, incx,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_dfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_csyr( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
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
	}

	bl1_csyr_blas( uplo,
	               m,
	               alpha,
	               x, incx,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_zsyr( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs )
{
	int       m_save    = m;
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
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
	}

	bl1_zsyr_blas( uplo,
	               m,
	               alpha,
	               x, incx,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}



// --- Classic routine wrappers ---

void bl1_ssyr_blas( uplo1_t uplo, int m, float* alpha, float* x, int incx, float* a, int lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssyr( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssyr( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

void bl1_dsyr_blas( uplo1_t uplo, int m, double* alpha, double* x, int incx, double* a, int lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsyr( cblas_order,
	            cblas_uplo,
	            m,
	            *alpha,
	            x, incx,
	            a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsyr( &blas_uplo,
	          &m,
	          alpha,
	          x, &incx,
	          a, &lda );
#endif
}

void bl1_csyr_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int lda )
{
	scomplex* x_copy;
	scomplex  beta;
	int       k   = 1;
	int       ldx = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &cblas_trans );

	x_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_csyrk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             alpha,
	             x_copy, ldx,
	             &beta,
	             a,      lda );

	bl1_cfree( x_copy );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &blas_trans );

	x_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_csyrk ( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            x_copy, &ldx,
	            &beta,
	            a,      &lda );

	bl1_cfree( x_copy );
#endif
}

void bl1_zsyr_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int lda )
{
	dcomplex* x_copy;
	dcomplex  beta;
	int       k   = 1;
	int       ldx = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &cblas_trans );

	x_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_zsyrk( cblas_order,
	             cblas_uplo,
	             cblas_trans,
	             m,
	             k,
	             alpha,
	             x_copy, ldx,
	             &beta,
	             a,      lda );

	bl1_zfree( x_copy );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &blas_trans );

	x_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_zsyrk ( &blas_uplo,
	            &blas_trans,
	            &m,
	            &k,
	            alpha,
	            x_copy, &ldx,
	            &beta,
	            a,      &lda );

	bl1_zfree( x_copy );
#endif
}

