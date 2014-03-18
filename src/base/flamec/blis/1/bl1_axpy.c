
#include "blis1.h"

void bl1_saxpy( int n, float* alpha, float* x, int incx, float* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_saxpy( n,
	             *alpha,
	             x, incx,
	             y, incy );
#else
	F77_saxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_daxpy( int n, double* alpha, double* x, int incx, double* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_daxpy( n,
	             *alpha,
	             x, incx,
	             y, incy );
#else
	F77_daxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_caxpy( int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_caxpy( n,
	             alpha,
	             x, incx,
	             y, incy );
#else
	F77_caxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

void bl1_zaxpy( int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zaxpy( n,
	             alpha,
	             x, incx,
	             y, incy );
#else
	F77_zaxpy( &n,
	           alpha,
	           x, &incx,
	           y, &incy );
#endif
}

