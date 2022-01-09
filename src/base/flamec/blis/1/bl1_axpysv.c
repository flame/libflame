/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_saxpysv( integer n, float* alpha0, float* alpha1, float* x, integer incx, float* beta, float* y, integer incy )
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

void bl1_daxpysv( integer n, double* alpha0, double* alpha1, double* x, integer incx, double* beta, double* y, integer incy )
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

void bl1_caxpysv( integer n, scomplex* alpha0, scomplex* alpha1, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy )
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

void bl1_zaxpysv( integer n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy )
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

