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

/*
   Effective computation:

     rho_xz = beta * rho_xz + x * z;
     rho_yz = beta * rho_yz + y * z;
     rho_wz = beta * rho_wz + w * z;

   where x, y, and w are optionally conjugated.
*/

void bl1_sdotsv3( conj1_t    conjxyw,
                  integer       n,
                  float*    x, integer inc_x,
                  float*    y, integer inc_y,
                  float*    w, integer inc_w,
                  float*    z, integer inc_z,
                  float*    beta,
                  float*    rho_xz,
                  float*    rho_yz,
                  float*    rho_wz )
{
	bl1_abort();
}


void bl1_ddotsv3( conj1_t    conjxyw,
                  integer       n,
                  double*   x, integer inc_x,
                  double*   y, integer inc_y,
                  double*   w, integer inc_w,
                  double*   z, integer inc_z,
                  double*   beta,
                  double*   rho_xz,
                  double*   rho_yz,
                  double*   rho_wz )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict x1;
	double*   restrict y1;
	double*   restrict w1;
	double*   restrict z1;
	double             rho1, rho2, rho3;
	double             x1c, y1c, w1c, z1c;
	integer                i;

	integer                n_pre;
	integer                n_run;
	integer                n_left;

	v2df_t             rho1v, rho2v, rho3v;
	v2df_t             x1v, y1v, w1v, z1v;
	v2df_t             x2v, y2v, w2v, z2v;
	
	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_w != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	if ( ( unsigned long ) z % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) y % 16 == 0 ||
		     ( unsigned long ) w % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	x1 = x;
	y1 = y;
	w1 = w;
	z1 = z;

	rho1 = 0.0;
	rho2 = 0.0;
	rho3 = 0.0;

	if ( n_pre == 1 )
	{
		x1c = *x1;
		y1c = *y1;
		w1c = *w1;
		z1c = *z1;

		rho1 += x1c * z1c;
		rho2 += y1c * z1c;
		rho3 += w1c * z1c;

		x1 += inc_x;
		y1 += inc_y;
		w1 += inc_w;
		z1 += inc_z;
	}

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();
	rho3v.v = _mm_setzero_pd();

	for ( i = 0; i < n_run; ++i )
	{
		x1v.v = _mm_load_pd( ( double* )x1 );
		y1v.v = _mm_load_pd( ( double* )y1 );
		w1v.v = _mm_load_pd( ( double* )w1 );
		z1v.v = _mm_load_pd( ( double* )z1 );

		rho1v.v += x1v.v * z1v.v;
		rho2v.v += y1v.v * z1v.v;
		rho3v.v += w1v.v * z1v.v;

		x2v.v = _mm_load_pd( ( double* )(x1 + 2) );
		y2v.v = _mm_load_pd( ( double* )(y1 + 2) );
		w2v.v = _mm_load_pd( ( double* )(w1 + 2) );
		z2v.v = _mm_load_pd( ( double* )(z1 + 2) );

		rho1v.v += x2v.v * z2v.v;
		rho2v.v += y2v.v * z2v.v;
		rho3v.v += w2v.v * z2v.v;

		x1 += 4;
		y1 += 4;
		w1 += 4;
		z1 += 4;
	}

	rho1 += rho1v.d[0] + rho1v.d[1];
	rho2 += rho2v.d[0] + rho2v.d[1];
	rho3 += rho3v.d[0] + rho3v.d[1];

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			w1c = *w1;
			z1c = *z1;

			rho1 += x1c * z1c;
			rho2 += y1c * z1c;
			rho3 += w1c * z1c;

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}
	}

	*rho_xz = *beta * *rho_xz + rho1;
	*rho_yz = *beta * *rho_yz + rho2;
	*rho_wz = *beta * *rho_wz + rho3;
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict x1;
	double*   restrict y1;
	double*   restrict w1;
	double*   restrict z1;
	double             rho1, rho2, rho3;
	double             x1c, y1c, w1c, z1c;
	double             x2c, y2c, w2c, z2c;
	integer                i;

	integer                n_pre;
	integer                n_run;
	integer                n_left;
	
	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_w != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	//if ( ( unsigned long ) z % 16 != 0 )
	//{
	//	if ( ( unsigned long ) x % 16 == 0 ||
	//	     ( unsigned long ) w % 16 == 0 ||
	//	     ( unsigned long ) y % 16 == 0 ) bl1_abort();
	//
	//	n_pre = 1;
	//}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	x1 = x;
	y1 = y;
	w1 = w;
	z1 = z;

	rho1 = 0.0;
	rho2 = 0.0;
	rho3 = 0.0;

	if ( n_pre == 1 )
	{
		x1c = *x1;
		y1c = *y1;
		w1c = *w1;
		z1c = *z1;

		rho1 += x1c * z1c;
		rho2 += y1c * z1c;
		rho3 += w1c * z1c;

		x1 += inc_x;
		y1 += inc_y;
		w1 += inc_w;
		z1 += inc_z;
	}

	for ( i = 0; i < n_run; ++i )
	{
		x1c = *x1;
		x2c = *(x1 + 1);
		y1c = *y1;
		y2c = *(y1 + 1);
		w1c = *w1;
		w2c = *(w1 + 1);
		z1c = *z1;
		z2c = *(z1 + 1);

		rho1 += x1c * z1c + x2c * z2c;
		rho2 += y1c * z1c + y2c * z2c;
		rho3 += w1c * z1c + w2c * z2c;

		x1 += 2*inc_x;
		y1 += 2*inc_y;
		w1 += 2*inc_w;
		z1 += 2*inc_z;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			w1c = *w1;
			z1c = *z1;

			rho1 += x1c * z1c;
			rho2 += y1c * z1c;
			rho3 += w1c * z1c;

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}
	}

	*rho_xz = *beta * *rho_xz + rho1;
	*rho_yz = *beta * *rho_yz + rho2;
	*rho_wz = *beta * *rho_wz + rho3;
}
#endif


void bl1_cdotsv3( conj1_t    conjxyw,
                  integer       n,
                  scomplex* x, integer inc_x,
                  scomplex* y, integer inc_y,
                  scomplex* w, integer inc_w,
                  scomplex* z, integer inc_z,
                  scomplex* beta,
                  scomplex* rho_xz,
                  scomplex* rho_yz,
                  scomplex* rho_wz )
{
	bl1_abort();
}


void bl1_zdotsv3( conj1_t    conjxyw,
                  integer       n,
                  dcomplex* x, integer inc_x,
                  dcomplex* y, integer inc_y,
                  dcomplex* w, integer inc_w,
                  dcomplex* z, integer inc_z,
                  dcomplex* beta,
                  dcomplex* rho_xz,
                  dcomplex* rho_yz,
                  dcomplex* rho_wz )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict x1;
	dcomplex* restrict y1;
	dcomplex* restrict w1;
	dcomplex* restrict z1;
	integer                i;
	v2df_t r1v, rho1v;
	v2df_t r2v, rho2v;
	v2df_t r3v, rho3v;
	v2df_t z11v, z12v;
	v2df_t x1v, x1rv;
	v2df_t y1v, y1rv;
	v2df_t w1v, w1rv;

	x1 = x;
	y1 = y;
	w1 = w;
	z1 = z;

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();
	rho3v.v = _mm_setzero_pd();

	if ( bl1_is_conj( conjxyw ) )
	{
		v2df_t bcac, adbd;

		for ( i = 0; i < n; ++i )
		{
			z11v.v = _mm_loaddup_pd( ( double* )&(z1->real) );
			z12v.v = _mm_loaddup_pd( ( double* )&(z1->imag) );

			x1v.v  = _mm_load_pd( ( double* )x1 );
			x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
			bcac.v = x1rv.v * z11v.v;
			adbd.v = x1v.v  * z12v.v;
			rho1v.v = rho1v.v + _mm_addsub_pd( bcac.v, adbd.v );

			y1v.v  = _mm_load_pd( ( double* )y1 );
			y1rv.v = _mm_shuffle_pd( y1v.v, y1v.v, _MM_SHUFFLE2 (0,1) );
			bcac.v = y1rv.v * z11v.v;
			adbd.v = y1v.v  * z12v.v;
			rho2v.v = rho2v.v + _mm_addsub_pd( bcac.v, adbd.v );

			w1v.v  = _mm_load_pd( ( double* )w1 );
			w1rv.v = _mm_shuffle_pd( w1v.v, w1v.v, _MM_SHUFFLE2 (0,1) );
			bcac.v = w1rv.v * z11v.v;
			adbd.v = w1v.v  * z12v.v;
			rho3v.v = rho3v.v + _mm_addsub_pd( bcac.v, adbd.v );

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}

		rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );
		rho2v.v = _mm_shuffle_pd( rho2v.v, rho2v.v, _MM_SHUFFLE2 (0,1) );
		rho3v.v = _mm_shuffle_pd( rho3v.v, rho3v.v, _MM_SHUFFLE2 (0,1) );

		rho1v.d[1] = -rho1v.d[1];
		rho2v.d[1] = -rho2v.d[1];
		rho3v.d[1] = -rho3v.d[1];
	}
	else
	{
		v2df_t cada, dbcb;

		for ( i = 0; i < n; ++i )
		{
			z11v.v = _mm_loaddup_pd( ( double* )&(z1->real) );
			z12v.v = _mm_loaddup_pd( ( double* )&(z1->imag) );

			x1v.v  = _mm_load_pd( ( double* )x1 );
			x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
			cada.v = x1v.v  * z11v.v;
			dbcb.v = x1rv.v * z12v.v;
			rho1v.v = rho1v.v + _mm_addsub_pd( cada.v, dbcb.v );

			y1v.v  = _mm_load_pd( ( double* )y1 );
			y1rv.v = _mm_shuffle_pd( y1v.v, y1v.v, _MM_SHUFFLE2 (0,1) );
			cada.v = y1v.v  * z11v.v;
			dbcb.v = y1rv.v * z12v.v;
			rho2v.v = rho2v.v + _mm_addsub_pd( cada.v, dbcb.v );

			w1v.v  = _mm_load_pd( ( double* )w1 );
			w1rv.v = _mm_shuffle_pd( w1v.v, w1v.v, _MM_SHUFFLE2 (0,1) );
			cada.v = w1v.v  * z11v.v;
			dbcb.v = w1rv.v * z12v.v;
			rho3v.v = rho3v.v + _mm_addsub_pd( cada.v, dbcb.v );

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}
	}

    //bl1_zscals( beta, rho_xz );
    //bl1_zscals( beta, rho_yz );
    //bl1_zscals( beta, rho_wz );
	{
		v2df_t ab, ba, cc, dd, acbc, bdad;

		ab.v = _mm_load_pd( ( double* )beta );
		ba.v = _mm_shuffle_pd( ab.v, ab.v, _MM_SHUFFLE2 (0,1) );

		cc.v = _mm_loaddup_pd( ( double* )&(rho_xz->real) );
		dd.v = _mm_loaddup_pd( ( double* )&(rho_xz->imag) );
		acbc.v = ab.v * cc.v;
		bdad.v = ba.v * dd.v;
		r1v.v = _mm_addsub_pd( acbc.v, bdad.v );

		cc.v = _mm_loaddup_pd( ( double* )&(rho_yz->real) );
		dd.v = _mm_loaddup_pd( ( double* )&(rho_yz->imag) );
		acbc.v = ab.v * cc.v;
		bdad.v = ba.v * dd.v;
		r2v.v = _mm_addsub_pd( acbc.v, bdad.v );

		cc.v = _mm_loaddup_pd( ( double* )&(rho_wz->real) );
		dd.v = _mm_loaddup_pd( ( double* )&(rho_wz->imag) );
		acbc.v = ab.v * cc.v;
		bdad.v = ba.v * dd.v;
		r3v.v = _mm_addsub_pd( acbc.v, bdad.v );
	}

	//rho_xz->real = rho_xz->real + rho1.real;
	//rho_xz->imag = rho_xz->imag + rho1.imag;
	rho1v.v = r1v.v + rho1v.v;
	_mm_store_pd( ( double* )rho_xz, rho1v.v );

	//rho_yz->real = rho_yz->real + rho2.real;
	//rho_yz->imag = rho_yz->imag + rho2.imag;
	rho2v.v = r2v.v + rho2v.v;
	_mm_store_pd( ( double* )rho_yz, rho2v.v );

	//rho_wz->real = rho_wz->real + rho3.real;
	//rho_wz->imag = rho_wz->imag + rho3.imag;
	rho3v.v = r3v.v + rho3v.v;
	_mm_store_pd( ( double* )rho_wz, rho3v.v );
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict x1;
	dcomplex* restrict y1;
	dcomplex* restrict w1;
	dcomplex* restrict z1;
	dcomplex           rho1, rho2, rho3;
	dcomplex           x1c, y1c, w1c, z1c;
	integer                i;

	x1 = x;
	y1 = y;
	w1 = w;
	z1 = z;

	rho1.real = 0.0; rho1.imag = 0.0;
	rho2.real = 0.0; rho2.imag = 0.0;
	rho3.real = 0.0; rho3.imag = 0.0;

	if ( bl1_is_conj( conjxyw ) )
	{
		for ( i = 0; i < n; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			w1c = *w1;
			z1c = *z1;

			rho1.real += x1c.real * z1c.real - -x1c.imag * z1c.imag;
			rho1.imag += x1c.real * z1c.imag + -x1c.imag * z1c.real;

			rho2.real += y1c.real * z1c.real - -y1c.imag * z1c.imag;
			rho2.imag += y1c.real * z1c.imag + -y1c.imag * z1c.real;

			rho3.real += w1c.real * z1c.real - -w1c.imag * z1c.imag;
			rho3.imag += w1c.real * z1c.imag + -w1c.imag * z1c.real;

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			w1c = *w1;
			z1c = *z1;

			rho1.real += x1c.real * z1c.real - x1c.imag * z1c.imag;
			rho1.imag += x1c.real * z1c.imag + x1c.imag * z1c.real;

			rho2.real += y1c.real * z1c.real - y1c.imag * z1c.imag;
			rho2.imag += y1c.real * z1c.imag + y1c.imag * z1c.real;

			rho3.real += w1c.real * z1c.real - w1c.imag * z1c.imag;
			rho3.imag += w1c.real * z1c.imag + w1c.imag * z1c.real;

			x1 += inc_x;
			y1 += inc_y;
			w1 += inc_w;
			z1 += inc_z;
		}
	}

    bl1_zscals( beta, rho_xz );
    bl1_zscals( beta, rho_yz );
    bl1_zscals( beta, rho_wz );

	rho_xz->real = rho_xz->real + rho1.real;
	rho_xz->imag = rho_xz->imag + rho1.imag;

	rho_yz->real = rho_yz->real + rho2.real;
	rho_yz->imag = rho_yz->imag + rho2.imag;

	rho_wz->real = rho_wz->real + rho3.real;
	rho_wz->imag = rho_wz->imag + rho3.imag;
}
#endif

