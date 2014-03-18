
#include "blis1.h"

void bl1_shemv( uplo1_t uplo, conj1_t conj, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	bl1_ssymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bl1_dhemv( uplo1_t uplo, conj1_t conj, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	bl1_dsymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bl1_chemv( uplo1_t uplo, conj1_t conj, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex  zero = bl1_c0();
	scomplex  one  = bl1_c1();
	scomplex* x_conj;
	scomplex* ax;
	int       lda, inca;
	int       incx_conj;
	int       incax;

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

	// We want to handle the case where A is conjugated, but without
	// explicitly or conjugating A. To do so, we leverage the fact that
	// computing the product conj(A) * x is equivalent to computing
	// conj( A * conj(x) ).
	if ( bl1_is_conj( conj ) )
	{
		// We need a temporary vector so we can create a conjugated copy of x.
		x_conj    = bl1_callocv( m );
		incx_conj = 1;

		bl1_ccopyv( BLIS1_CONJUGATE,
		            m,
		            x,      incx,
                    x_conj, incx_conj );

		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y (and x).
		ax    = bl1_callocv( m );
		incax = 1;
		
		// Compute A * conj(x) where x is the temporary copy of x created above.
		bl1_chemv_blas( uplo,
		                m,
                        &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax,     incax );

		// Scale y by beta.
		bl1_cscalv( BLIS1_NO_CONJUGATE,
                    m,
                    beta,
                    y, incy );

		// And finally, accumulate alpha * conj( A * conj(x) ) into y.
		bl1_caxpyv( BLIS1_CONJUGATE,
                    m,
		            alpha,
                    ax, incax,
                    y,  incy);

		// Free the temporary vectors for x and Ax.
		bl1_cfree( x_conj );
		bl1_cfree( ax );
	}
	else // noconj
	{
		bl1_chemv_blas( uplo,
		                m,
		                alpha,
		                a, lda,
		                x, incx,
		                beta,
		                y, incy );
	}

	// Free the temporary contiguous matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_zhemv( uplo1_t uplo, conj1_t conj, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex  zero = bl1_z0();
	dcomplex  one  = bl1_z1();
	dcomplex* x_conj;
	dcomplex* ax;
	int       lda, inca;
	int       incx_conj;
	int       incax;

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

	// We want to handle the case where A is conjugated, but without
	// explicitly or conjugating A. To do so, we leverage the fact that
	// computing the product conj(A) * x is equivalent to computing
	// conj( A * conj(x) ).
	if ( bl1_is_conj( conj ) )
	{
		// We need a temporary vector so we can create a conjugated copy of x.
		x_conj    = bl1_zallocv( m );
		incx_conj = 1;

		bl1_zcopyv( BLIS1_CONJUGATE,
		            m,
		            x,      incx,
                    x_conj, incx_conj );

		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y (and x).
		ax    = bl1_zallocv( m );
		incax = 1;
		
		// Compute A * conj(x) where x is the temporary copy of x created above.
		bl1_zhemv_blas( uplo,
		                m,
		                &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax,     incax );

		// Scale y by beta.
		bl1_zscalv( BLIS1_NO_CONJUGATE,
                    m,
                    beta,
                    y, incy );

		// And finally, accumulate alpha * conj( A * conj(x) ) into y.
		bl1_zaxpyv( BLIS1_CONJUGATE,
                    m,
                    alpha,
                    ax, incax,
                    y,  incy);

		// Free the temporary vectors for x and Ax.
		bl1_zfree( x_conj );
		bl1_zfree( ax );
	}
	else // noconj
	{
		bl1_zhemv_blas( uplo,
		                m,
		                alpha,
		                a, lda,
		                x, incx,
		                beta,
		                y, incy );
	}

	// Free the temporary contiguous matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_chemv_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_chemv( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_chemv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_zhemv_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zhemv( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zhemv( &blas_uplo,
	           &m,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

