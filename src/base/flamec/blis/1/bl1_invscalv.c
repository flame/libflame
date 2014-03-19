/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sinvscalv( conj1_t conj, int n, float* alpha, float* x, int incx )
{
	float alpha_inv;

	if ( bl1_seq1( alpha ) ) return;

	alpha_inv = 1.0F / *alpha;

	bl1_sscal( n,
	           &alpha_inv,
	           x, incx );
}

void bl1_dinvscalv( conj1_t conj, int n, double* alpha, double* x, int incx )
{
	double alpha_inv;

	if ( bl1_deq1( alpha ) ) return;

	alpha_inv = 1.0 / *alpha;

	bl1_dscal( n,
	           &alpha_inv,
	           x, incx );
}

void bl1_csinvscalv( conj1_t conj, int n, float* alpha, scomplex* x, int incx )
{
	float alpha_inv;

	if ( bl1_seq1( alpha ) ) return;

	alpha_inv = 1.0F / *alpha;

	bl1_csscal( n,
	            &alpha_inv,
	            x, incx );
}

void bl1_cinvscalv( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx )
{
	scomplex alpha_inv;

	if ( bl1_ceq1( alpha ) ) return;

	bl1_cinvert2s( conj, alpha, &alpha_inv );

	bl1_cscal( n,
	           &alpha_inv,
	           x, incx );
}

void bl1_zdinvscalv( conj1_t conj, int n, double* alpha, dcomplex* x, int incx )
{
	double alpha_inv;

	if ( bl1_deq1( alpha ) ) return;

	alpha_inv = 1.0 / *alpha;

	bl1_zdscal( n,
	            &alpha_inv,
	            x, incx );
}

void bl1_zinvscalv( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx )
{
	dcomplex alpha_inv;

	if ( bl1_zeq1( alpha ) ) return;

	bl1_zinvert2s( conj, alpha, &alpha_inv );

	bl1_zscal( n,
	           &alpha_inv,
	           x, incx );
}

