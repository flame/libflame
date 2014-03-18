
#include "blis1.h"

void bl1_snrm2( int n, float* x, int incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_snrm2( n,
	                     x, incx );
#else
	*norm = F77_snrm2( &n,
	                   x, &incx );
#endif
}

void bl1_dnrm2( int n, double* x, int incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dnrm2( n,
	                     x, incx );
#else
	*norm = F77_dnrm2( &n,
	                   x, &incx );
#endif
}

void bl1_cnrm2( int n, scomplex* x, int incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_scnrm2( n,
	                      x, incx );
#else
	*norm = F77_scnrm2( &n,
	                    x, &incx );
#endif
}

void bl1_znrm2( int n, dcomplex* x, int incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dznrm2( n,
	                      x, incx );
#else
	*norm = F77_dznrm2( &n,
	                    x, &incx );
#endif
}

