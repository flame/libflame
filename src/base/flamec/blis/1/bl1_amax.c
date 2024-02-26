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

void bl1_samax( integer n, float* x, integer incx, integer* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_isamax( n,
	                       x, incx );
#else
	*index = F77_isamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_damax( integer n, double* x, integer incx, integer* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_idamax( n,
	                       x, incx );
#else
	*index = F77_idamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_camax( integer n, scomplex* x, integer incx, integer* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_icamax( n,
	                       x, incx );
#else
	*index = F77_icamax( &n,
	                     x, &incx ) - 1;
#endif
}

void bl1_zamax( integer n, dcomplex* x, integer incx, integer* index )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	*index = cblas_izamax( n,
	                       x, incx );
#else
	*index = F77_izamax( &n,
	                     x, &incx ) - 1;
#endif
}

