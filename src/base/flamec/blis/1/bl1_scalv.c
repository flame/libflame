
#include "blis1.h"

void bl1_sscalv( conj1_t conj, int n, float* alpha, float* x, int incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	bl1_sscal( n,
	           alpha,
	           x, incx );
}

void bl1_dscalv( conj1_t conj, int n, double* alpha, double* x, int incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	bl1_dscal( n,
	           alpha,
	           x, incx );
}

void bl1_csscalv( conj1_t conj, int n, float* alpha, scomplex* x, int incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_seq1( alpha ) ) return;

	bl1_csscal( n,
	            alpha,
	            x, incx );
}

void bl1_cscalv( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx )
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

void bl1_zdscalv( conj1_t conj, int n, double* alpha, dcomplex* x, int incx )
{
	// Return early if possible.
	if ( bl1_zero_dim1( n ) ) return;
	if ( bl1_deq1( alpha ) ) return;

	bl1_zdscal( n,
	            alpha,
	            x, incx );
}

void bl1_zscalv( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx )
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

