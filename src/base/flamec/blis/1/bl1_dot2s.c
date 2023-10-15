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

void bl1_sdot2s( conj1_t conj, integer n, float* alpha, float* x, integer incx, float* y, integer incy, float* beta, float* rho )
{
	float dot;

	bl1_sdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0F * (*alpha) * dot;
}

void bl1_ddot2s( conj1_t conj, integer n, double* alpha, double* x, integer incx, double* y, integer incy, double* beta, double* rho )
{
	double dot;

	bl1_ddot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0 * (*alpha) * dot;
}

void bl1_cdot2s( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* beta, scomplex* rho )
{
	scomplex dotxy;
	scomplex dotyx;
	scomplex alpha_d    = *alpha;
	scomplex alphac_d   = *alpha;
	scomplex beta_d     = *beta;
	scomplex rho_d      = *rho;

	alphac_d.imag *= -1.0F;

	bl1_cdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dotxy );

	bl1_cdot( conj,
	          n,
	          y, incy,
	          x, incx,
	          &dotyx );

	rho->real = beta_d.real   * rho_d.real - beta_d.imag   * rho_d.imag +
	            alpha_d.real  * dotxy.real - alpha_d.imag  * dotxy.imag +
	            alphac_d.real * dotyx.real - alphac_d.imag * dotyx.imag; 
	rho->imag = beta_d.real   * rho_d.imag + beta_d.imag   * rho_d.real +
	            alpha_d.real  * dotxy.imag + alpha_d.imag  * dotxy.real +
	            alphac_d.real * dotyx.imag + alphac_d.imag * dotyx.real; 
}

void bl1_zdot2s( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* beta, dcomplex* rho )
{
	dcomplex dotxy;
	dcomplex dotyx;
	dcomplex alpha_d    = *alpha;
	dcomplex alphac_d   = *alpha;
	dcomplex beta_d     = *beta;
	dcomplex rho_d      = *rho;

	alphac_d.imag *= -1.0;

	bl1_zdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dotxy );

	bl1_zdot( conj,
	          n,
	          y, incy,
	          x, incx,
	          &dotyx );

	rho->real = beta_d.real   * rho_d.real - beta_d.imag   * rho_d.imag +
	            alpha_d.real  * dotxy.real - alpha_d.imag  * dotxy.imag +
	            alphac_d.real * dotyx.real - alphac_d.imag * dotyx.imag; 
	rho->imag = beta_d.real   * rho_d.imag + beta_d.imag   * rho_d.real +
	            alpha_d.real  * dotxy.imag + alpha_d.imag  * dotxy.real +
	            alphac_d.real * dotyx.imag + alphac_d.imag * dotyx.real; 
}

