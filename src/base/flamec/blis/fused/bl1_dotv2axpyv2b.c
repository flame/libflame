/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

/*
   Effective computation:

     rho1 = conj(a1) * x;
     w    = w + kappa1 * a1;
     rho2 = conj(a2) * x;
     w    = w + kappa2 * a2;
*/

void bl1_sdotv2axpyv2b( integer       n,
                        float*    a1, integer inc_a1,
                        float*    a2, integer inc_a2,
                        float*    x,  integer inc_x,
                        float*    kappa1,
                        float*    kappa2,
                        float*    rho1,
                        float*    rho2,
                        float*    w, integer inc_w )
{
	bl1_abort();
}


void bl1_ddotv2axpyv2b( integer       n,
                        double*   a1, integer inc_a1,
                        double*   a2, integer inc_a2,
                        double*   x,  integer inc_x,
                        double*   kappa1,
                        double*   kappa2,
                        double*   rho1,
                        double*   rho2,
                        double*   w, integer inc_w )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict alpha1;
	double*   restrict alpha2;
	double*   restrict chi1;
	double*   restrict omega1;
	double             rho1_c;
	double             rho2_c;
	integer                i;

	integer                n_pre;
	integer                n_run;
	integer                n_left;
	
	v2df_t    k1v, rho1v;
	v2df_t    k2v, rho2v;
	v2df_t    a11v, a12v, x1v, w1v;
	v2df_t    a21v, a22v, x2v, w2v;
	
	if ( inc_a1 != 1 ||
	     inc_a2 != 1 ||
	     inc_x  != 1 ||
	     inc_w  != 1 ) bl1_abort();

	n_pre = 0;
	if ( ( unsigned long ) a1 % 16 != 0 )
	{
		if ( ( unsigned long ) a2 % 16 == 0 ||
		     ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) w % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	alpha1   = a1;
	alpha2   = a2;
	chi1     = x;
	omega1   = w;

	rho1_c = 0.0;
	rho2_c = 0.0;

	if ( n_pre == 1 )
	{
		double   kappa1_c = *kappa1;
		double   kappa2_c = *kappa2;
		double   alpha1_c   = *alpha1;
		double   alpha2_c   = *alpha2;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		rho1_c   += alpha1_c * chi1_c;
		omega1_c += kappa1_c * alpha1_c;

		rho2_c   += alpha2_c * chi1_c;
		omega1_c += kappa2_c * alpha2_c;

		*omega1 = omega1_c;

		alpha1   += inc_a1;
		alpha2   += inc_a2;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();

	k1v.v = _mm_loaddup_pd( ( double* )kappa1 );
	k2v.v = _mm_loaddup_pd( ( double* )kappa2 );

	for ( i = 0; i < n_run; ++i )
	{
		a11v.v = _mm_load_pd( ( double* )alpha1 );
		a12v.v = _mm_load_pd( ( double* )alpha2 );
		x1v.v  = _mm_load_pd( ( double* )chi1 );
		w1v.v  = _mm_load_pd( ( double* )omega1 );

		rho1v.v += a11v.v * x1v.v;
		w1v.v += k1v.v * a11v.v;

		rho2v.v += a12v.v * x1v.v;
		w1v.v += k2v.v * a12v.v;

		_mm_store_pd( ( double* )omega1, w1v.v );

		a21v.v = _mm_load_pd( ( double* )(alpha1 + 2) );
		a22v.v = _mm_load_pd( ( double* )(alpha2 + 2) );
		x2v.v  = _mm_load_pd( ( double* )(chi1 + 2) );
		w2v.v  = _mm_load_pd( ( double* )(omega1 + 2) );

		rho1v.v += a21v.v * x2v.v;
		w2v.v += k1v.v * a21v.v;

		rho2v.v += a22v.v * x2v.v;
		w2v.v += k2v.v * a22v.v;

		_mm_store_pd( ( double* )(omega1 + 2), w2v.v );

		alpha1   += 4;
		alpha2   += 4;
		chi1     += 4;
		omega1   += 4;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			double   kappa1_c = *kappa1;
			double   kappa2_c = *kappa2;
			double   alpha1_c   = *alpha1;
			double   alpha2_c   = *alpha2;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			rho1_c   += alpha1_c * chi1_c;
			omega1_c += kappa1_c * alpha1_c;

			rho2_c   += alpha2_c * chi1_c;
			omega1_c += kappa2_c * alpha2_c;

			*omega1 = omega1_c;

			alpha1   += inc_a1;
			alpha2   += inc_a2;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	rho1_c += rho1v.d[0] + rho1v.d[1];
	rho2_c += rho2v.d[0] + rho2v.d[1];

	*rho1 = rho1_c;
	*rho2 = rho2_c;
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict alpha1;
	double*   restrict alpha2;
	double*   restrict chi1;
	double*   restrict omega1;
	double             kappa1_c;
	double             kappa2_c;
	double             rho1_c;
	double             rho2_c;
	integer                i;

	integer                n_pre;
	integer                n_run;
	integer                n_left;
	
	if ( inc_a1 != 1 ||
	     inc_a2 != 1 ||
	     inc_x  != 1 ||
	     inc_w  != 1 ) bl1_abort();

	n_pre = 0;
	//if ( ( unsigned long ) a1 % 16 != 0 )
	//{
	//	if ( ( unsigned long ) a2 % 16 == 0 ||
	//	     ( unsigned long ) x % 16 == 0 ||
	//	     ( unsigned long ) w % 16 == 0 ) bl1_abort();
	//
	//	n_pre = 1;
	//}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	alpha1   = a1;
	alpha2   = a2;
	chi1     = x;
	omega1   = w;

	rho1_c = 0.0;
	rho2_c = 0.0;

	kappa1_c = *kappa1;
	kappa2_c = *kappa2;

	if ( n_pre == 1 )
	{
		double   alpha1_c   = *alpha1;
		double   alpha2_c   = *alpha2;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		rho1_c   += alpha1_c * chi1_c;
		omega1_c += kappa1_c * alpha1_c;

		rho2_c   += alpha2_c * chi1_c;
		omega1_c += kappa2_c * alpha2_c;

		*omega1 = omega1_c;

		alpha1   += inc_a1;
		alpha2   += inc_a2;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	for ( i = 0; i < n_run; ++i )
	{
		double   alpha11_c   = *alpha1;
		double   alpha21_c   = *(alpha1 + 1);
		double   alpha12_c   = *alpha2;
		double   alpha22_c   = *(alpha2 + 1);
		double   chi1_c     = *chi1;
		double   chi2_c     = *(chi1 + 1);
		double   omega1_c   = *omega1;
		double   omega2_c   = *(omega1 + 1);

		// rho1 += conj(alpha1) * chi1;
		rho1_c += alpha11_c * chi1_c;
		rho1_c += alpha21_c * chi2_c;

		// omega1 += kappa1 * alpha1;
		omega1_c += kappa1_c * alpha11_c;
		omega2_c += kappa1_c * alpha21_c;

		// rho2 += conj(alpha2) * chi1;
		rho2_c += alpha12_c * chi1_c;
		rho2_c += alpha22_c * chi2_c;

		// omega1 += kappa2 * alpha2;
		omega1_c += kappa2_c * alpha12_c;
		omega2_c += kappa2_c * alpha22_c;

		*omega1       = omega1_c;
		*(omega1 + 1) = omega2_c;

		alpha1   += 2*inc_a1;
		alpha2   += 2*inc_a2;
		chi1     += 2*inc_x;
		omega1   += 2*inc_w;
	}

	if ( n_left > 0 )
	{
		for ( i = 0; i < n_left; ++i )
		{
			double   alpha1_c   = *alpha1;
			double   alpha2_c   = *alpha2;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			rho1_c   += alpha1_c * chi1_c;
			omega1_c += kappa1_c * alpha1_c;

			rho2_c   += alpha2_c * chi1_c;
			omega1_c += kappa2_c * alpha2_c;

			*omega1 = omega1_c;

			alpha1   += inc_a1;
			alpha2   += inc_a2;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	*rho1 = rho1_c;
	*rho2 = rho2_c;
}
#endif


void bl1_cdotv2axpyv2b( integer       n,
                        scomplex* a1, integer inc_a1,
                        scomplex* a2, integer inc_a2,
                        scomplex* x,  integer inc_x,
                        scomplex* kappa1,
                        scomplex* kappa2,
                        scomplex* rho1,
                        scomplex* rho2,
                        scomplex* w, integer inc_w )
{
	bl1_abort();
}


void bl1_zdotv2axpyv2b( integer       n,
                        dcomplex* a1, integer inc_a1,
                        dcomplex* a2, integer inc_a2,
                        dcomplex* x,  integer inc_x,
                        dcomplex* kappa1,
                        dcomplex* kappa2,
                        dcomplex* rho1,
                        dcomplex* rho2,
                        dcomplex* w, integer inc_w )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict alpha1;
	dcomplex* restrict alpha2;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	integer                i;

	v2df_t    kappa1v, kappa1rv;
	v2df_t    kappa2v, kappa2rv;
	v2df_t    rho1v;
	v2df_t    rho2v;
	v2df_t    a11v, a12v;
	v2df_t    a21v, a22v;
	v2df_t    x1v, x1rv;
	v2df_t    w1v;
	v2df_t    acbc, bdad;
	v2df_t    adac, bcbd;

	if ( inc_a1 != 1 ||
	     inc_a2 != 1 ||
	     inc_x  != 1 ||
	     inc_w  != 1 ) bl1_abort();

	alpha1   = a1;
	alpha2   = a2;
	chi1     = x;
	omega1   = w;

	rho1v.v = _mm_setzero_pd();
	rho2v.v = _mm_setzero_pd();

	kappa1v.v  = _mm_load_pd( ( double* )kappa1 );
	kappa1rv.v = _mm_shuffle_pd( kappa1v.v, kappa1v.v, _MM_SHUFFLE2 (0,1) );
	kappa2v.v  = _mm_load_pd( ( double* )kappa2 );
	kappa2rv.v = _mm_shuffle_pd( kappa2v.v, kappa2v.v, _MM_SHUFFLE2 (0,1) );

	for ( i = 0; i < n; ++i )
	{
		//dcomplex omega1_c = *omega1;
		w1v.v  = _mm_load_pd( ( double* )omega1 );

		//dcomplex chi1_c   = *chi1;
		x1v.v  = _mm_load_pd( ( double* )chi1 );


		//dcomplex alpha1_c = *alpha1;
		a11v.v  = _mm_loaddup_pd( ( double* )&(alpha1->real) );
		a12v.v  = _mm_loaddup_pd( ( double* )&(alpha1->imag) );

		//rho1_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		//rho1_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;
        x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
        adac.v = a11v.v * x1rv.v;
        bcbd.v = a12v.v * x1v.v;
        rho1v.v = rho1v.v + _mm_addsub_pd( adac.v, bcbd.v );

		//omega1_c.real += kappa1_c.real * alpha1_c.real - kappa1_c.imag * alpha1_c.imag;
		//omega1_c.imag += kappa1_c.real * alpha1_c.imag + kappa1_c.imag * alpha1_c.real;
		acbc.v = kappa1v.v  * a11v.v;
		bdad.v = kappa1rv.v * a12v.v;
		w1v.v += _mm_addsub_pd( acbc.v, bdad.v );


		//dcomplex alpha2_c = *alpha2;
		a21v.v  = _mm_loaddup_pd( ( double* )&(alpha2->real) );
		a22v.v  = _mm_loaddup_pd( ( double* )&(alpha2->imag) );

		//rho2_c.real += alpha2_c.real * chi1_c.real - -alpha2_c.imag * chi1_c.imag;
		//rho2_c.imag += alpha2_c.real * chi1_c.imag + -alpha2_c.imag * chi1_c.real;
        x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
        adac.v = a21v.v * x1rv.v;
        bcbd.v = a22v.v * x1v.v;
        rho2v.v = rho2v.v + _mm_addsub_pd( adac.v, bcbd.v );

		//omega1_c.real += kappa2_c.real * alpha2_c.real - kappa2_c.imag * alpha2_c.imag;
		//omega1_c.imag += kappa2_c.real * alpha2_c.imag + kappa2_c.imag * alpha2_c.real;
		acbc.v = kappa2v.v  * a21v.v;
		bdad.v = kappa2rv.v * a22v.v;
		w1v.v += _mm_addsub_pd( acbc.v, bdad.v );


		//*omega1 = omega1_c;
		_mm_store_pd( ( double* )omega1, w1v.v );


		//alpha1   += inc_a1;
		//alpha2   += inc_a2;
		//chi1     += inc_x;
		//omega1   += inc_w;
		alpha1   += 1;
		alpha2   += 1;
		chi1     += 1;
		omega1   += 1;
	}

	rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );
	rho2v.v = _mm_shuffle_pd( rho2v.v, rho2v.v, _MM_SHUFFLE2 (0,1) );

	//*rho1 = rho1_c;
	//*rho2 = rho2_c;
	_mm_store_pd( ( double* )rho1, rho1v.v );
	_mm_store_pd( ( double* )rho2, rho2v.v );
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict alpha1;
	dcomplex* restrict alpha2;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	dcomplex           rho1_c;
	dcomplex           rho2_c;
	dcomplex           kappa1_c;
	dcomplex           kappa2_c;
	integer                i;

	alpha1   = a1;
	alpha2   = a2;
	chi1     = x;
	omega1   = w;

	rho1_c.real = 0.0;
	rho1_c.imag = 0.0;
	rho2_c.real = 0.0;
	rho2_c.imag = 0.0;

	kappa1_c = *kappa1;
	kappa1_c = *kappa1;
	kappa2_c = *kappa2;
	kappa2_c = *kappa2;

	for ( i = 0; i < n; ++i )
	{
		dcomplex alpha1_c   = *alpha1;
		dcomplex alpha2_c   = *alpha2;
		dcomplex chi1_c     = *chi1;
		dcomplex omega1_c   = *omega1;

		// rho1 += conj(alpha1) * chi1;
		rho1_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		rho1_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;

		// omega1 += kappa1 * alpha1;
		omega1_c.real += kappa1_c.real * alpha1_c.real - kappa1_c.imag * alpha1_c.imag;
		omega1_c.imag += kappa1_c.real * alpha1_c.imag + kappa1_c.imag * alpha1_c.real;

		// rho2 += conj(alpha2) * chi1;
		rho2_c.real += alpha2_c.real * chi1_c.real - -alpha2_c.imag * chi1_c.imag;
		rho2_c.imag += alpha2_c.real * chi1_c.imag + -alpha2_c.imag * chi1_c.real;

		// omega1 += kappa2 * alpha2;
		omega1_c.real += kappa2_c.real * alpha2_c.real - kappa2_c.imag * alpha2_c.imag;
		omega1_c.imag += kappa2_c.real * alpha2_c.imag + kappa2_c.imag * alpha2_c.real;

		*omega1 = omega1_c;

		alpha1   += inc_a1;
		alpha2   += inc_a2;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	*rho1 = rho1_c;
	*rho2 = rho2_c;
}
#endif

