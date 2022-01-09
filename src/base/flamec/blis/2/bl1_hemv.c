/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_shemv( uplo1_t uplo, conj1_t conj, integer m, float* alpha, float* a, integer a_rs, integer a_cs, float* x, integer incx, float* beta, float* y, integer incy )
{
	bl1_ssymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bl1_dhemv( uplo1_t uplo, conj1_t conj, integer m, double* alpha, double* a, integer a_rs, integer a_cs, double* x, integer incx, double* beta, double* y, integer incy )
{
	bl1_dsymv( uplo,
	           m,
	           alpha,
	           a, a_rs, a_cs,
	           x, incx,
	           beta,
	           y, incy );
}

void bl1_chemv( uplo1_t uplo, conj1_t conj, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
{
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	scomplex  zero = bl1_c0();
	scomplex  one  = bl1_c1();
	scomplex* x_conj;
	scomplex* ax;
	integer       lda, inca;
	integer       incx_conj;
	integer       incax;

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

void bl1_zhemv( uplo1_t uplo, conj1_t conj, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
{
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	dcomplex  zero = bl1_z0();
	dcomplex  one  = bl1_z1();
	dcomplex* x_conj;
	dcomplex* ax;
	integer       lda, inca;
	integer       incx_conj;
	integer       incax;

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

void bl1_chemv_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* a, integer lda, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
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

void bl1_zhemv_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
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

