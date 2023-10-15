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

void bl1_sdot( conj1_t conj, integer n, float* x, integer incx, float* y, integer incy, float* rho )
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

void bl1_ddot( conj1_t conj, integer n, double* x, integer incx, double* y, integer incy, double* rho )
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

void bl1_cdot( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* rho )
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

void bl1_zdot( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* rho )
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

void bl1_cdot_in( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* rho )
{
	scomplex* xip;
	scomplex* yip;
	scomplex  xi;
	scomplex  yi;
	scomplex  rho_temp;
	integer       i;

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

void bl1_zdot_in( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* rho )
{
	dcomplex* xip;
	dcomplex* yip;
	dcomplex  xi;
	dcomplex  yi;
	dcomplex  rho_temp;
	integer       i;

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

