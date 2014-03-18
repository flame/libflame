
#include "blis1.h"

void bl1_sdot2s( conj1_t conj, int n, float* alpha, float* x, int incx, float* y, int incy, float* beta, float* rho )
{
	float dot;

	bl1_sdot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0F * (*alpha) * dot;
}

void bl1_ddot2s( conj1_t conj, int n, double* alpha, double* x, int incx, double* y, int incy, double* beta, double* rho )
{
	double dot;

	bl1_ddot( conj,
	          n,
	          x, incx,
	          y, incy,
	          &dot );

	*rho = (*beta) * (*rho) + 2.0 * (*alpha) * dot;
}

void bl1_cdot2s( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho )
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

void bl1_zdot2s( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho )
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

