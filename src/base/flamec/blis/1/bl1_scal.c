
#include "blis1.h"

void bl1_sscal( int n, float* alpha, float* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_sscal( n,
	             *alpha,
	             x, incx );
#else
	F77_sscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_dscal( int n, double* alpha, double* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_dscal( n,
	             *alpha,
	             x, incx );
#else
	F77_dscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_csscal( int n, float* alpha, scomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_csscal( n,
	              *alpha,
	              x, incx );
#else
	F77_csscal( &n,
	            alpha,
	            x, &incx );
#endif
}

void bl1_cscal( int n, scomplex* alpha, scomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_cscal( n,
	             alpha,
	             x, incx );
#else
	F77_cscal( &n,
	           alpha,
	           x, &incx );
#endif
}

void bl1_zdscal( int n, double* alpha, dcomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zdscal( n,
	              *alpha,
	              x, incx );
#else
	F77_zdscal( &n,
	            alpha,
	            x, &incx );
#endif
}

void bl1_zscal( int n, dcomplex* alpha, dcomplex* x, int incx )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zscal( n,
	             alpha,
	             x, incx );
#else
	F77_zscal( &n,
	           alpha,
	           x, &incx );
#endif
}

