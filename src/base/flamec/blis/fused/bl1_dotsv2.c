
#include "blis1.h"

/*
   Effective computation:

     rho_xz = beta * rho_xz + x * z;
     rho_yz = beta * rho_yz + y * z;

   where x and y are optionally conjugated.
*/

void bl1_sdotsv2( conj1_t    conjxy,
                  int       n,
                  float*    x, int inc_x,
                  float*    y, int inc_y,
                  float*    z, int inc_z,
                  float*    beta,
                  float*    rho_xz,
                  float*    rho_yz )
{
	bl1_abort();
}


void bl1_ddotsv2( conj1_t    conjxy,
                  int       n,
                  double*   x, int inc_x,
                  double*   y, int inc_y,
                  double*   z, int inc_z,
                  double*   beta,
                  double*   rho_xz,
                  double*   rho_yz )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict x1;
	double*   restrict y1;
	double*   restrict z1;
	double             rho1, rho2;
	double             x1c, y1c, z1c;
	int                i;

	int                n_pre;
	int                n_run;
	int                n_left;

	v2df_t             rho1v, rho2v;
	v2df_t             x1v, y1v, z1v;
	v2df_t             x2v, y2v, z2v;
	
	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	if ( ( unsigned long ) z % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) y % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	x1 = x;
	y1 = y;
	z1 = z;

	rho1 = 0.0;
	rho2 = 0.0;

	if ( n_pre == 1 )
	{
		x1c = *x1;
		y1c = *y1;
		z1c = *z1;

		rho1 += x1c * z1c;
		rho2 += y1c * z1c;

		x1 += inc_x;
		y1 += inc_y;
		z1 += inc_z;
	}

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();

	for ( i = 0; i < n_run; ++i )
	{
		x1v.v = _mm_load_pd( ( double* )x1 );
		y1v.v = _mm_load_pd( ( double* )y1 );
		z1v.v = _mm_load_pd( ( double* )z1 );

		x2v.v = _mm_load_pd( ( double* )(x1 + 2) );
		y2v.v = _mm_load_pd( ( double* )(y1 + 2) );
		z2v.v = _mm_load_pd( ( double* )(z1 + 2) );

		rho1v.v += x1v.v * z1v.v;
		rho2v.v += y1v.v * z1v.v;

		rho1v.v += x2v.v * z2v.v;
		rho2v.v += y2v.v * z2v.v;

		x1 += 4;
		y1 += 4;
		z1 += 4;
	}

	rho1 += rho1v.d[0] + rho1v.d[1];
	rho2 += rho2v.d[0] + rho2v.d[1];

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			z1c = *z1;

			rho1 += x1c * z1c;
			rho2 += y1c * z1c;

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}
	}

	*rho_xz = *beta * *rho_xz + rho1;
	*rho_yz = *beta * *rho_yz + rho2;
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict x1;
	double*   restrict y1;
	double*   restrict z1;
	double             rho1, rho2;
	double             x1c, y1c, z1c;
	double             x2c, y2c, z2c;
	int                i;

	int                n_pre;
	int                n_run;
	int                n_left;
	
	if ( inc_x != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	//if ( ( unsigned long ) z % 16 != 0 )
	//{
	//	if ( ( unsigned long ) x % 16 == 0 ||
	//	     ( unsigned long ) y % 16 == 0 ) bl1_abort();
	//
	//	n_pre = 1;
	//}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	x1 = x;
	y1 = y;
	z1 = z;

	rho1 = 0.0;
	rho2 = 0.0;

	if ( n_pre == 1 )
	{
		x1c = *x1;
		y1c = *y1;
		z1c = *z1;

		rho1 += x1c * z1c;
		rho2 += y1c * z1c;

		x1 += inc_x;
		y1 += inc_y;
		z1 += inc_z;
	}

	for ( i = 0; i < n_run; ++i )
	{
		x1c = *x1;
		x2c = *(x1 + 1);
		y1c = *y1;
		y2c = *(y1 + 1);
		z1c = *z1;
		z2c = *(z1 + 1);

		rho1 += x1c * z1c + x2c * z2c;
		rho2 += y1c * z1c + y2c * z2c;

		x1 += 2*inc_x;
		y1 += 2*inc_y;
		z1 += 2*inc_z;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			z1c = *z1;

			rho1 += x1c * z1c;
			rho2 += y1c * z1c;

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}
	}

	*rho_xz = *beta * *rho_xz + rho1;
	*rho_yz = *beta * *rho_yz + rho2;
}
#endif


void bl1_cdotsv2( conj1_t    conjxy,
                  int       n,
                  scomplex* x, int inc_x,
                  scomplex* y, int inc_y,
                  scomplex* z, int inc_z,
                  scomplex* beta,
                  scomplex* rho_xz,
                  scomplex* rho_yz )
{
	bl1_abort();
}


void bl1_zdotsv2( conj1_t    conjxy,
                  int       n,
                  dcomplex* x, int inc_x,
                  dcomplex* y, int inc_y,
                  dcomplex* z, int inc_z,
                  dcomplex* beta,
                  dcomplex* rho_xz,
                  dcomplex* rho_yz )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict x1;
	dcomplex* restrict y1;
	dcomplex* restrict z1;
	int                i;
	v2df_t r1v, rho1v;
	v2df_t r2v, rho2v;
	v2df_t z11v, z12v;
	v2df_t x1v, x1rv;
	v2df_t y1v, y1rv;

	x1 = x;
	y1 = y;
	z1 = z;

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();

	if ( bl1_is_conj( conjxy ) )
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

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}

		rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );
		rho2v.v = _mm_shuffle_pd( rho2v.v, rho2v.v, _MM_SHUFFLE2 (0,1) );

		rho1v.d[1] = -rho1v.d[1];
		rho2v.d[1] = -rho2v.d[1];
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

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}
	}

    //bl1_zscals( beta, rho_xz );
    //bl1_zscals( beta, rho_yz );
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
	}

	//rho_xz->real = rho_xz->real + rho1.real;
	//rho_xz->imag = rho_xz->imag + rho1.imag;
	rho1v.v = r1v.v + rho1v.v;
	_mm_store_pd( ( double* )rho_xz, rho1v.v );

	//rho_yz->real = rho_yz->real + rho2.real;
	//rho_yz->imag = rho_yz->imag + rho2.imag;
	rho2v.v = r2v.v + rho2v.v;
	_mm_store_pd( ( double* )rho_yz, rho2v.v );
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict x1;
	dcomplex* restrict y1;
	dcomplex* restrict z1;
	dcomplex           rho1, rho2;
	dcomplex           x1c, y1c, z1c;
	int                i;

	x1 = x;
	y1 = y;
	z1 = z;

	rho1.real = 0.0; rho1.imag = 0.0;
	rho2.real = 0.0; rho2.imag = 0.0;

	if ( bl1_is_conj( conjxy ) )
	{
		for ( i = 0; i < n; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			z1c = *z1;

			rho1.real += x1c.real * z1c.real - -x1c.imag * z1c.imag;
			rho1.imag += x1c.real * z1c.imag + -x1c.imag * z1c.real;

			rho2.real += y1c.real * z1c.real - -y1c.imag * z1c.imag;
			rho2.imag += y1c.real * z1c.imag + -y1c.imag * z1c.real;

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			x1c = *x1;
			y1c = *y1;
			z1c = *z1;

			rho1.real += x1c.real * z1c.real - x1c.imag * z1c.imag;
			rho1.imag += x1c.real * z1c.imag + x1c.imag * z1c.real;

			rho2.real += y1c.real * z1c.real - y1c.imag * z1c.imag;
			rho2.imag += y1c.real * z1c.imag + y1c.imag * z1c.real;

			x1 += inc_x;
			y1 += inc_y;
			z1 += inc_z;
		}
	}

    bl1_zscals( beta, rho_xz );
    bl1_zscals( beta, rho_yz );

	rho_xz->real = rho_xz->real + rho1.real;
	rho_xz->imag = rho_xz->imag + rho1.imag;

	rho_yz->real = rho_yz->real + rho2.real;
	rho_yz->imag = rho_yz->imag + rho2.imag;
}
#endif

