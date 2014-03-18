
#include "blis1.h"

void bl1_samax( int n, float* x, int incx, int* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_isamax( n,
	                       x, incx );
#else
	*index = F77_isamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_damax( int n, double* x, int incx, int* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_idamax( n,
	                       x, incx );
#else
	*index = F77_idamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_camax( int n, scomplex* x, int incx, int* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_icamax( n,
	                       x, incx );
#else
	*index = F77_icamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_zamax( int n, dcomplex* x, int incx, int* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_izamax( n,
	                       x, incx );
#else
	*index = F77_izamax( &n,
	                     x, &incx ) - 1;
#endif
}

