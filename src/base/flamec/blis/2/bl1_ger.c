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

void bl1_sger( conj1_t conjx, conj1_t conjy, integer m, integer n, float* alpha, float* x, integer incx, float* y, integer incy, float* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	float*    a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		bl1_swap_ints( incx, incy );
		bl1_swap_conj( conjx, conjy );
		bl1_sswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	bl1_sger_blas( m,
	               n,
	               alpha,
	               x, incx,
	               y, incy,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_sfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_dger( conj1_t conjx, conj1_t conjy, integer m, integer n, double* alpha, double* x, integer incx, double* y, integer incy, double* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	double*   a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		bl1_swap_ints( incx, incy );
		bl1_swap_conj( conjx, conjy );
		bl1_dswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	bl1_dger_blas( m,
	               n,
	               alpha,
	               x, incx,
	               y, incy,
	               a, lda );

	// Free the temporary contiguous matrix.
	bl1_dfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_cger( conj1_t conjx, conj1_t conjy, integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	scomplex* x_conj;
	integer       incx_conj;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		bl1_swap_ints( incx, incy );
		bl1_swap_conj( conjx, conjy );
		bl1_cswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated.
	if ( bl1_is_conj( conjx ) )
	{
		x_conj    = bl1_callocv( m );
		incx_conj = 1;

		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	// Conjugation of y is supported in the BLAS.
	if ( bl1_is_conj( conjy ) )
	{
		bl1_cgerc_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}
	else
	{
		bl1_cgeru_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conjx ) )
		bl1_cfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_cfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_zger( conj1_t conjx, conj1_t conjy, integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	integer       n_save    = n;
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	dcomplex* x_conj;
	integer       incx_conj;
	integer       lda, inca;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

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
		bl1_swap_ints( incx, incy );
		bl1_swap_conj( conjx, conjy );
		bl1_zswap_pointers( x, y );
	}

	// Initialize with values assuming no conjugation of x.
	x_conj    = x;
	incx_conj = incx;

	// We need a temporary vector for the cases when x is conjugated.
	if ( bl1_is_conj( conjx ) )
	{
		x_conj    = bl1_zallocv( m );
		incx_conj = 1;

		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );
	}

	// Conjugation of y is supported in the BLAS.
	if ( bl1_is_conj( conjy ) )
	{
		bl1_zgerc_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}
	else
	{
		bl1_zgeru_blas( m,
		                n,
		                alpha,
		                x_conj, incx_conj,
		                y,      incy, 
		                a,      lda );
	}

	// Free the temporary conjugated x vector.
	if ( bl1_is_conj( conjx ) )
		bl1_zfree( x_conj );

	// Free the temporary contiguous matrix.
	bl1_zfree_saved_contigm( m_save,
	                         n_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_sger_blas( integer m, integer n, float* alpha, float* x, integer incx, float* y, integer incy, float* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_sger( cblas_order,
	            m,
	            n,
	            *alpha,
	            x, incx,
	            y, incy,
	            a, lda );
#else
	F77_sger( &m,
	          &n,
	          alpha,
	          x, &incx,
	          y, &incy,
	          a, &lda );
#endif
}

void bl1_dger_blas( integer m, integer n, double* alpha, double* x, integer incx, double* y, integer incy, double* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_dger( cblas_order,
	            m,
	            n,
	            *alpha,
	            x, incx,
	            y, incy,
	            a, lda );
#else
	F77_dger( &m,
	          &n,
	          alpha,
	          x, &incx,
	          y, &incy,
	          a, &lda );
#endif
}

void bl1_cgerc_blas( integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_cgerc( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_cgerc ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bl1_cgeru_blas( integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_cgeru( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_cgeru ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bl1_zgerc_blas( integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_zgerc( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_zgerc ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

void bl1_zgeru_blas( integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;

	cblas_zgeru( cblas_order,
	             m,
	             n,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	F77_zgeru ( &m,
	            &n,
	            alpha,
	            x, &incx,
	            y, &incy,
	            a, &lda );
#endif
}

