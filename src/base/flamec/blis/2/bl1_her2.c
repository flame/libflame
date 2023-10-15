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

void bl1_sher2( uplo1_t uplo, conj1_t conj, integer m, float* alpha, float* x, integer incx, float* y, integer incy, float* a, integer a_rs, integer a_cs )
{
	bl1_ssyr2( uplo,
	           m,
	           alpha,
	           x, incx,
	           y, incy,
	           a, a_rs, a_cs );
}

void bl1_dher2( uplo1_t uplo, conj1_t conj, integer m, double* alpha, double* x, integer incx, double* y, integer incy, double* a, integer a_rs, integer a_cs )
{
	bl1_dsyr2( uplo,
	           m,
	           alpha,
	           x, incx,
	           y, incy,
	           a, a_rs, a_cs );
}

void bl1_cher2( uplo1_t uplo, conj1_t conj, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	scomplex* x_conj;
	scomplex* y_conj;
	integer       incx_conj;
	integer       incy_conj;
	integer       lda, inca;

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

	// Initialize with values assuming no conjugation of ( x * y' ) or
	// ( y * x' ).
	x_conj    = x;
	incx_conj = incx;
	y_conj    = y;
	incy_conj = incy;

	// We want to handle the case where ( x * y' ) and ( y * x' ) are
	// conjugated, but without explicitly conjugating the matrices. To do
	// so, we leverage the fact that computing the products conj( x * y' )
	// and conj( y * x' ) is equivalent to computing ( conj(x) * conj(y)' )
	// and ( conj(y) * conj(x)' ), respectively.
	if ( bl1_is_conj( conj ) )
	{
		x_conj    = bl1_callocv( m );
		incx_conj = 1;

		y_conj    = bl1_callocv( m );
		incy_conj = 1;

		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );

		bl1_ccopyv( BLIS1_CONJUGATE,
                    m,
                    y,      incy,
                    y_conj, incy_conj );
	}

	bl1_cher2_blas( uplo,
	                m,
	                alpha,
	                x_conj, incx_conj,
	                y_conj, incy_conj,
	                a,      lda );

	// Free the temporary conjugated x and y vectors.
	if ( bl1_is_conj( conj ) )
	{
		bl1_cfree( x_conj );
		bl1_cfree( y_conj );
	}

	// Free the temporary contiguous matrix.
	bl1_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_zher2( uplo1_t uplo, conj1_t conj, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	dcomplex* x_conj;
	dcomplex* y_conj;
	integer       incx_conj;
	integer       incy_conj;
	integer       lda, inca;

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

	// Initialize with values assuming no conjugation of ( x * y' ) or
	// ( y * x' ).
	x_conj    = x;
	incx_conj = incx;
	y_conj    = y;
	incy_conj = incy;

	// We want to handle the case where ( x * y' ) and ( y * x' ) are
	// conjugated, but without explicitly conjugating the matrices. To do
	// so, we leverage the fact that computing the products conj( x * y' )
	// and conj( y * x' ) is equivalent to computing ( conj(x) * conj(y)' )
	// and ( conj(y) * conj(x)' ), respectively.
	if ( bl1_is_conj( conj ) )
	{
		x_conj    = bl1_zallocv( m );
		incx_conj = 1;

		y_conj    = bl1_zallocv( m );
		incy_conj = 1;

		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    x,      incx,
                    x_conj, incx_conj );

		bl1_zcopyv( BLIS1_CONJUGATE,
                    m,
                    y,      incy,
                    y_conj, incy_conj );
	}

	bl1_zher2_blas( uplo,
	                m,
	                alpha,
	                x_conj, incx_conj,
	                y_conj, incy_conj,
	                a,      lda );

	// Free the temporary conjugated x and y vectors.
	if ( bl1_is_conj( conj ) )
	{
		bl1_zfree( x_conj );
		bl1_zfree( y_conj );
	}

	// Free the temporary contiguous matrix.
	bl1_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_cher2_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_cher2( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_cher2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bl1_zher2_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zher2( cblas_order,
	             cblas_uplo,
	             m,
	             alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zher2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

