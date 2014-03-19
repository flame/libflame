/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

/*
   Effective computation:

     rho_xz = beta * rho_xz + x * z;
     rho_yz = beta * rho_yz + y * z;

   where x and y are optionally conjugated.
*/

dcomplex zdotc_( int*      n,
                 dcomplex* x, int* inc_x,
                 dcomplex* z, int* inc_z )
{
	dcomplex* restrict x1;
	dcomplex* restrict z1;
	int                i;
	v2df_t rho1v;
	v2df_t z11v, z12v;
	v2df_t x1v, x1rv;
	dcomplex rho;
	int    n1 = *n;
	int    incx = *inc_x;
	int    incz = *inc_z;

	x1 = x;
	z1 = z;

	rho1v.v = _mm_setzero_pd();

	{
		v2df_t bcac, adbd;

		for ( i = 0; i < n1; ++i )
		{
			z11v.v = _mm_loaddup_pd( ( double* )&(z1->real) );
			z12v.v = _mm_loaddup_pd( ( double* )&(z1->imag) );

			x1v.v  = _mm_load_pd( ( double* )x1 );
			x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
			bcac.v = x1rv.v * z11v.v;
			adbd.v = x1v.v  * z12v.v;
			rho1v.v = rho1v.v + _mm_addsub_pd( bcac.v, adbd.v );

			x1 += incx;
			z1 += incz;
		}

		rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );

		rho1v.d[1] = -rho1v.d[1];
	}

	rho.real = rho1v.d[0];
	rho.imag = rho1v.d[1];

	return rho;
}

