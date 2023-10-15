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

void bl1_sinvert2s( conj1_t conj, float* alpha, float* beta )
{
	float  one = 1.0F;

	*beta = one / *alpha;
}

void bl1_dinvert2s( conj1_t conj, double* alpha, double* beta )
{
	double one = 1.0;

	*beta = one / *alpha;
}

void bl1_cinvert2s( conj1_t conj, scomplex* alpha, scomplex* beta )
{
	float  temp;
	float  s, xr_s, xi_s;

	s           = bl1_fmaxabs( alpha->real, alpha->imag ); \
	xr_s        = alpha->real / s;
	xi_s        = alpha->imag / s;
	temp        = xr_s * alpha->real + xi_s * alpha->imag;

	beta->real =  xr_s / temp;
	beta->imag = -xi_s / temp;

	if ( bl1_is_conj( conj ) )
		bl1_cconjs( beta );
}

void bl1_zinvert2s( conj1_t conj, dcomplex* alpha, dcomplex* beta )
{
	double temp;
	double s, xr_s, xi_s;

	s           = bl1_fmaxabs( alpha->real, alpha->imag ); \
	xr_s        = alpha->real / s;
	xi_s        = alpha->imag / s;
	temp        = xr_s * alpha->real + xi_s * alpha->imag;

	beta->real =  xr_s / temp;
	beta->imag = -xi_s / temp;

	if ( bl1_is_conj( conj ) )
		bl1_zconjs( beta );
}

