
#include "blis1.h"

void bl1_scopy( int m, float* x, int incx, float* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_scopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_scopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_dcopy( int m, double* x, int incx, double* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_dcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_dcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_ccopy( int m, scomplex* x, int incx, scomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_ccopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_ccopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

void bl1_zcopy( int m, dcomplex* x, int incx, dcomplex* y, int incy )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	cblas_zcopy( m,
	             x, incx, 
	             y, incy );
#else
	F77_zcopy( &m,
	           x, &incx, 
	           y, &incy );
#endif
}

