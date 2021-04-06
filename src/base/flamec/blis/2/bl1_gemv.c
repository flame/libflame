/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sgemv( trans1_t transa, conj1_t conjx, integer m, integer n, float* alpha, float* a, integer a_rs, integer a_cs, float* x, integer incx, float* beta, float* y, integer incy )
{
	float*    a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) )
	{
		integer n_elem;

		if ( bl1_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bl1_sscalv( BLIS1_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_screate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( m, n );
		bl1_swap_ints( lda, inca );
		bl1_toggle_trans( transa );
	}

	bl1_sgemv_blas( transa,
	                m,
	                n,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_dgemv( trans1_t transa, conj1_t conjx, integer m, integer n, double* alpha, double* a, integer a_rs, integer a_cs, double* x, integer incx, double* beta, double* y, integer incy )
{
	double*   a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) )
	{
		integer n_elem;

		if ( bl1_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bl1_dscalv( BLIS1_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_dcreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( m, n );
		bl1_swap_ints( lda, inca );
		bl1_toggle_trans( transa );
	}

	bl1_dgemv_blas( transa,
	                m,
	                n,
	                alpha,
	                a, lda,
	                x, incx,
	                beta,
	                y, incy );

	// Free the temporary contiguous matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_cgemv( trans1_t transa, conj1_t conjx, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
{
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	scomplex  zero = bl1_c0();
	scomplex  one  = bl1_c1();
	scomplex* x_conj;
	scomplex* ax;
	integer       lda, inca;
	integer       n_x;
	integer       incx_conj;
	integer       incax;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) )
	{
		integer n_elem;

		if ( bl1_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bl1_cscalv( BLIS1_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_ccreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( m, n );
		bl1_swap_ints( lda, inca );
		bl1_toggle_trans( transa );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated, and
	// also for the cases where A is conjugated.
	if ( bl1_is_conj( conjx ) || bl1_is_conjnotrans( transa ) )
	{
		if ( bl1_does_trans( transa ) ) n_x = m;
		else                            n_x = n;

		x_conj    = bl1_callocv( n_x );
		incx_conj = 1;

		bl1_ccopyv( conjx,
		            n_x,
		            x,      incx,
                    x_conj, incx_conj );
	}

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	if ( bl1_is_conjnotrans( transa ) )
	{
		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y. We know we are not transposing, so y is length m.
		ax    = bl1_callocv( m );
		incax = 1;
		
		// Start by conjugating the contents of the temporary copy of x.
		bl1_cconjv( n,
		            x_conj, incx_conj );

		// Compute A * conj(x) where x is the temporary copy of x created above.
		bl1_cgemv_blas( BLIS1_NO_TRANSPOSE,
		                m,
		                n,
		                &one,
		                a,      lda,
		                x_conj, incx_conj,
		                &zero,
		                ax, incax );

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

		// Free the temporary vector for Ax.
		bl1_cfree( ax );
	}
	else // notrans, trans, or conjtrans
	{
		bl1_cgemv_blas( transa,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                x_conj, incx_conj,
		                beta,
		                y, incy );
	}

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conjx ) || bl1_is_conjnotrans( transa ) )
		bl1_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_zgemv( trans1_t transa, conj1_t conjx, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
{
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	dcomplex  zero = bl1_z0();
	dcomplex  one  = bl1_z1();
	dcomplex* x_conj;
	dcomplex* ax;
	integer       lda, inca;
	integer       n_x;
	integer       incx_conj;
	integer       incax;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) )
	{
		integer n_elem;

		if ( bl1_does_trans( transa ) ) n_elem = n;
		else                            n_elem = m;

		bl1_zscalv( BLIS1_NO_CONJUGATE,
		            n_elem,
		            beta,
		            y, incy );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_zcreate_contigm( m,
	                     n,
	                     a_save, a_rs_save, a_cs_save,
	                     &a,     &a_rs,     &a_cs );

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;

	// If A is a row-major matrix, then we can use the underlying column-major
	// BLAS implementation by fiddling with the parameters.
	if ( bl1_is_row_storage( a_rs, a_cs ) )
	{
		bl1_swap_ints( m, n );
		bl1_swap_ints( lda, inca );
		bl1_toggle_trans( transa );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated, and
	// also for the cases where A is conjugated.
	if ( bl1_is_conj( conjx ) || bl1_is_conjnotrans( transa ) )
	{
		if ( bl1_does_trans( transa ) ) n_x = m;
		else                            n_x = n;

		x_conj    = bl1_zallocv( n_x );
		incx_conj = 1;

		bl1_zcopyv( conjx,
		            n_x,
		            x,      incx,
                    x_conj, incx_conj );
	}

	// We want to handle the conjnotrans case, but without explicitly
	// conjugating A. To do so, we leverage the fact that computing the
	// product conj(A) * x is equivalent to computing conj( A * conj(x) ).
	if ( bl1_is_conjnotrans( transa ) )
	{
		// We need a temporary vector for the product A * conj(x), which is
		// conformal to y. We know we are not transposing, so y is length m.
		ax    = bl1_zallocv( m );
		incax = 1;
		
		// Start by conjugating the contents of the temporary copy of x.
		bl1_zconjv( n,
		            x_conj, incx_conj );

		// Compute A * conj(x) where x is the temporary copy of x created above.
		bl1_zgemv_blas( BLIS1_NO_TRANSPOSE,
		                m,
		                n,
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

		// Free the temporary vector for Ax.
		bl1_zfree( ax );
	}
	else // notrans, trans, or conjtrans
	{
		bl1_zgemv_blas( transa,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                x_conj, incx_conj,
		                beta,
		                y,      incy );
	}

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conjx ) || bl1_is_conjnotrans( transa ) )
		bl1_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_sgemv_blas( trans1_t transa, integer m, integer n, float* alpha, float* a, integer lda, float* x, integer incx, float* beta, float* y, integer incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_sgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_transa;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );

	F77_sgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_dgemv_blas( trans1_t transa, integer m, integer n, double* alpha, double* a, integer lda, double* x, integer incx, double* beta, double* y, integer incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_dgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             x, incx,
	             *beta,
	             y, incy );
#else
	char blas_transa;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );

	F77_dgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_cgemv_blas( trans1_t transa, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_cgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_transa;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );

	F77_cgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

void bl1_zgemv_blas( trans1_t transa, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );

	cblas_zgemv( cblas_order,
	             cblas_transa,
	             m,
	             n,
	             alpha,
	             a, lda,
	             x, incx,
	             beta,
	             y, incy );
#else
	char blas_transa;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );

	F77_zgemv( &blas_transa,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           x, &incx,
	           beta,
	           y, &incy );
#endif
}

