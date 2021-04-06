/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sscalv( conj1_t conj, integer n, float* alpha, float* x, integer incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	bl1_sscal( n,
	           alpha,
	           x, incx );
}

void bl1_dscalv( conj1_t conj, integer n, double* alpha, double* x, integer incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	bl1_dscal( n,
	           alpha,
	           x, incx );
}

void bl1_csscalv( conj1_t conj, integer n, float* alpha, scomplex* x, integer incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	bl1_csscal( n,
	            alpha,
	            x, incx );
}

void bl1_cscalv( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx )
{
	scomplex alpha_conj;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_ceq1( alpha ) ) return;

	bl1_ccopys( conj, alpha, &alpha_conj );

	bl1_cscal( n,
	           &alpha_conj,
	           x, incx );
}

void bl1_zdscalv( conj1_t conj, integer n, double* alpha, dcomplex* x, integer incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	bl1_zdscal( n,
	            alpha,
	            x, incx );
}

void bl1_zscalv( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx )
{
	dcomplex alpha_conj;

	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_zeq1( alpha ) ) return;

	bl1_zcopys( conj, alpha, &alpha_conj );

	bl1_zscal( n,
	           &alpha_conj,
	           x, incx );
}

