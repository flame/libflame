/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sdots( conj1_t conj, int n, float* alpha, float* x, int incx, float* y, int incy, float* beta, float* rho )
{
	float dot_prod;

	bl1_sdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	*rho = (*beta) * (*rho) + (*alpha) * dot_prod;
}

void bl1_ddots( conj1_t conj, int n, double* alpha, double* x, int incx, double* y, int incy, double* beta, double* rho )
{
	double dot_prod;

	bl1_ddot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	*rho = (*beta) * (*rho) + (*alpha) * dot_prod;
}

void bl1_cdots( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho )
{
	scomplex rho_orig = *rho;
	scomplex dot_prod;

	bl1_cdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	rho->real = beta->real  * rho_orig.real - beta->imag  * rho_orig.imag +
	            alpha->real * dot_prod.real - alpha->imag * dot_prod.imag;
	rho->imag = beta->real  * rho_orig.imag + beta->imag  * rho_orig.real +
	            alpha->real * dot_prod.imag + alpha->imag * dot_prod.real;
}

void bl1_zdots( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho )
{
	dcomplex rho_orig = *rho;
	dcomplex dot_prod;

	bl1_zdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot_prod );

	rho->real = beta->real  * rho_orig.real - beta->imag  * rho_orig.imag +
	            alpha->real * dot_prod.real - alpha->imag * dot_prod.imag;
	rho->imag = beta->real  * rho_orig.imag + beta->imag  * rho_orig.real +
	            alpha->real * dot_prod.imag + alpha->imag * dot_prod.real;
}

