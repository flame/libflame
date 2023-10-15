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

void bl1_ssyr2( uplo1_t uplo, integer m, float* alpha, float* x, integer incx, float* y, integer incy, float* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	float*    a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

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

	bl1_ssyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bl1_sfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_dsyr2( uplo1_t uplo, integer m, double* alpha, double* x, integer incx, double* y, integer incy, double* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	double*   a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	integer       lda, inca;

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

	bl1_dsyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bl1_dfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_csyr2( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
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
	}

	bl1_csyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bl1_cfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

void bl1_zsyr2( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs )
{
	integer       m_save    = m;
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
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
	}

	bl1_zsyr2_blas( uplo,
	                m,
	                alpha,
	                x, incx,
	                y, incy,
	                a, lda );

	// Free the temporary contiguous matrix.
	bl1_zfree_saved_contigm( m_save,
	                         m_save,
	                         a_save, a_rs_save, a_cs_save,
	                         &a,     &a_rs,     &a_cs );
}

// --- Classic routine wrappers ---

void bl1_ssyr2_blas( uplo1_t uplo, integer m, float* alpha, float* x, integer incx, float* y, integer incy, float* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssyr2( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssyr2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bl1_dsyr2_blas( uplo1_t uplo, integer m, double* alpha, double* x, integer incx, double* y, integer incy, double* a, integer lda )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsyr2( cblas_order,
	             cblas_uplo,
	             m,
	             *alpha,
	             x, incx,
	             y, incy,
	             a, lda );
#else
	char blas_uplo;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsyr2( &blas_uplo,
	           &m,
	           alpha,
	           x, &incx,
	           y, &incy,
	           a, &lda );
#endif
}

void bl1_csyr2_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda )
{
	scomplex* x_copy;
	scomplex* y_copy;
	scomplex  beta;
	integer       k   = 1;
	integer       ldx = m;
	integer       ldy = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &cblas_trans );

	x_copy = bl1_callocv( m );
	y_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_csyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              x_copy, ldx,
	              y_copy, ldy,
	              &beta,
	              a,      lda );

	bl1_cfree( x_copy );
	bl1_cfree( y_copy );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &blas_trans );

	x_copy = bl1_callocv( m );
	y_copy = bl1_callocv( m );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_csyr2k ( &blas_uplo,
	             &blas_trans,
	             &m,
	             &k,
	             alpha,
	             x_copy, &ldx,
	             y_copy, &ldy,
	             &beta,
	             a,      &lda );

	bl1_cfree( x_copy );
	bl1_cfree( y_copy );
#endif
}

void bl1_zsyr2_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda )
{
	dcomplex* x_copy;
	dcomplex* y_copy;
	dcomplex  beta;
	integer       k   = 1;
	integer       ldx = m;
	integer       ldy = m;

#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_UPLO      cblas_uplo;
	enum CBLAS_TRANSPOSE cblas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &cblas_trans );

	x_copy = bl1_zallocv( m );
	y_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	cblas_zsyr2k( cblas_order,
	              cblas_uplo,
	              cblas_trans,
	              m,
	              k,
	              alpha,
	              x_copy, ldx,
	              y_copy, ldy,
	              &beta,
	              a,      lda );

	bl1_zfree( x_copy );
	bl1_zfree( y_copy );
#else
	char blas_uplo;
	char blas_trans;

	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );
	bl1_param_map_to_netlib_trans( BLIS1_NO_TRANSPOSE, &blas_trans );

	x_copy = bl1_zallocv( m );
	y_copy = bl1_zallocv( m );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_copy, 1 );

	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            y,      incy,
	            y_copy, 1 );

	beta.real = 1.0;
	beta.imag = 0.0;

	F77_zsyr2k ( &blas_uplo,
	             &blas_trans,
	             &m,
	             &k,
	             alpha,
	             x_copy, &ldx,
	             y_copy, &ldy,
	             &beta,
	             a,      &lda );

	bl1_zfree( x_copy );
	bl1_zfree( y_copy );
#endif
}

