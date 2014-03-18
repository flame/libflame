
#include "blis1.h"

void bl1_sdot( conj1_t conj, int n, float* x, int incx, float* y, int incy, float* rho )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*rho = cblas_sdot( n,
	                   x, incx,
	                   y, incy );
#else
	*rho = F77_sdot( &n,
	                 x, &incx,
	                        y, &incy );
#endif
}

void bl1_ddot( conj1_t conj, int n, double* x, int incx, double* y, int incy, double* rho )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*rho = cblas_ddot( n,
	                   x, incx,
	                   y, incy );
#else
	*rho = F77_ddot( &n,
	                 x, &incx,
	                        y, &incy );
#endif
}

void bl1_cdot( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	if ( bl1_is_conj( conj ) )
	{
	    cblas_cdotc_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
	else // if ( !bl1_is_conj( conj ) )
	{
	    cblas_cdotu_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
#else
	bl1_cdot_in( conj,
	             n,
	             x, incx,
	             y, incy,
	             rho );
#endif
}

void bl1_zdot( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	if ( bl1_is_conj( conj ) )
	{
	    cblas_zdotc_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
	else // if ( !bl1_is_conj( conj ) )
	{
	    cblas_zdotu_sub( n,
		                 x, incx,
		                 y, incy,
		                 rho );
	}
#else
	bl1_zdot_in( conj,
	             n,
	             x, incx,
	             y, incy,
	             rho );
#endif
}


// --- Inlined helper implementations ---

void bl1_cdot_in( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho )
{
	scomplex* xip;
	scomplex* yip;
	scomplex  xi;
	scomplex  yi;
	scomplex  rho_temp;
	int       i;

	rho_temp.real = 0.0F;
	rho_temp.imag = 0.0F;
		
	xip = x;
	yip = y;
		
	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - -xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + -xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	else // if ( !bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	
	rho->real = rho_temp.real;
	rho->imag = rho_temp.imag;
}

void bl1_zdot_in( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho )
{
	dcomplex* xip;
	dcomplex* yip;
	dcomplex  xi;
	dcomplex  yi;
	dcomplex  rho_temp;
	int       i;

	rho_temp.real = 0.0;
	rho_temp.imag = 0.0;
		
	xip = x;
	yip = y;
		
	if ( bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - -xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + -xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	else // if ( !bl1_is_conj( conj ) )
	{
		for ( i = 0; i < n; ++i )
		{
			xi.real = xip->real;
			xi.imag = xip->imag;
			yi.real = yip->real;
			yi.imag = yip->imag;
			
			rho_temp.real += xi.real * yi.real - xi.imag * yi.imag;
			rho_temp.imag += xi.real * yi.imag + xi.imag * yi.real;

			xip += incx;
			yip += incy;
		}
	}
	
	rho->real = rho_temp.real;
	rho->imag = rho_temp.imag;
}

