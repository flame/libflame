
#include "blis1.h"

void bl1_saxpysv( int n, float* alpha0, float* alpha1, float* x, int incx, float* beta, float* y, int incy )
{
	float    alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bl1_sscal( n,
	           beta,
	           y, incy );

	bl1_saxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_daxpysv( int n, double* alpha0, double* alpha1, double* x, int incx, double* beta, double* y, int incy )
{
	double   alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod = (*alpha0) * (*alpha1);

	bl1_dscal( n,
	           beta,
	           y, incy );

	bl1_daxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_caxpysv( int n, scomplex* alpha0, scomplex* alpha1, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy )
{
	scomplex alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bl1_cscal( n,
	           beta,
	           y, incy );

	bl1_caxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

void bl1_zaxpysv( int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy )
{
	dcomplex alpha_prod;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;

	alpha_prod.real = alpha0->real * alpha1->real - alpha0->imag * alpha1->imag;
	alpha_prod.imag = alpha0->real * alpha1->imag + alpha0->imag * alpha1->real;

	bl1_zscal( n,
	           beta,
	           y, incy );

	bl1_zaxpy( n,
	           &alpha_prod,
	           x, incx,
	           y, incy );
}

