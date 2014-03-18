
#include "blis1.h"

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy )
{
	float*    a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	float*    x_temp;
    int       incx_temp;

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

void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy )
{
	double*   a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	double*   x_temp;
    int       incx_temp;

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

void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	scomplex* x_temp;
    int       incx_temp;

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

void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex* a_save    = a;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	dcomplex* x_temp;
    int       incx_temp;

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

