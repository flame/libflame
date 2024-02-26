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

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float* alpha, float* a, integer a_rs, integer a_cs, float* x, integer incx, float* beta, float* y, integer incy )
{
	float*    a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	float*    x_temp;
    integer       incx_temp;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_screate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bl1_sallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bl1_scopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bl1_strsv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bl1_sscalv( BLIS1_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bl1_saxpyv( BLIS1_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bl1_sfree( x_temp );

	// Free the temporary contiguous matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double* alpha, double* a, integer a_rs, integer a_cs, double* x, integer incx, double* beta, double* y, integer incy )
{
	double*   a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	double*   x_temp;
    integer       incx_temp;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_dcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bl1_dallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bl1_dcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bl1_dtrsv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bl1_dscalv( BLIS1_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bl1_daxpyv( BLIS1_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bl1_dfree( x_temp );

	// Free the temporary contiguous matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
{
	scomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	scomplex* x_temp;
    integer       incx_temp;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_ccreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bl1_callocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bl1_ccopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bl1_ctrsv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bl1_cscalv( BLIS1_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bl1_caxpyv( BLIS1_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bl1_cfree( x_temp );

	// Free the temporary contiguous matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
{
	dcomplex* a_save    = a;
	integer       a_rs_save = a_rs;
	integer       a_cs_save = a_cs;
	dcomplex* x_temp;
    integer       incx_temp;

	// Return early if possible.
	if ( bl1_zero_dim1( m ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of the matrix rather than the original matrix.
	bl1_zcreate_contigmr( uplo,
	                      m,
	                      m,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	// Allocate a temporary vector conformal to x.
	x_temp    = bl1_zallocv( m );
	incx_temp = 1;

	// Copy x to a temporary vector.
	bl1_zcopyv( BLIS1_NO_CONJUGATE,
	            m,
	            x,      incx,
	            x_temp, incx_temp );

	// Perform the operation, storing the result to x_temp.
	bl1_ztrsv( uplo,
	           trans,
	           diag,
	           m,
	           a,      a_rs, a_cs,
	           x_temp, incx_temp );

	// Scale y by beta.
	bl1_zscalv( BLIS1_NO_CONJUGATE,
	            m,
	            beta,
	            y, incy );

	// Axpy the partial result in x_temp into y.
	bl1_zaxpyv( BLIS1_NO_CONJUGATE,
	            m,
	            alpha,
	            x_temp, incx_temp,
	            y,      incy );

	// Free the temporary vector.
	bl1_zfree( x_temp );

	// Free the temporary contiguous matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );
}

