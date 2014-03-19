/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_sasum( int n, float* x, int incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_sasum( n,
	                     x, incx );
#else
	*norm = F77_sasum( &n,
	                   x, &incx );
#endif
}

void bl1_dasum( int n, double* x, int incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dasum( n,
	                     x, incx );
#else
	*norm = F77_dasum( &n,
	                   x, &incx );
#endif
}

void bl1_casum( int n, scomplex* x, int incx, float* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_scasum( n,
	                      x, incx );
#else
	*norm = F77_scasum( &n,
	                    x, &incx );
#endif
}

void bl1_zasum( int n, dcomplex* x, int incx, double* norm )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*norm = cblas_dzasum( n,
	                      x, incx );
#else
	*norm = F77_dzasum( &n,
	                    x, &incx );
#endif
}

