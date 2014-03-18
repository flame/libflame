
#include "FLAME.h"

/*
   Effective computation:

   y = y + alpha1 * x1;
   y = y + alpha2 * x2;
*/

void bl1_saxpyv2b( int       n,
                   float*    alpha1,
                   float*    alpha2,
                   float*    x1, int inc_x1,
                   float*    x2, int inc_x2,
                   float*    y,  int inc_y )
{
	bl1_abort();
}


void bl1_daxpyv2b( int       n,
                   double*   alpha1,
                   double*   alpha2,
                   double*   x1, int inc_x1,
                   double*   x2, int inc_x2,
                   double*   y,  int inc_y )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict chi2;
	double*   restrict psi1;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	v2df_t    a1v, a2v;
	v2df_t    x11v, x12v;
	v2df_t    x21v, x22v;
	v2df_t    y1v;
	v2df_t    y2v;

	if ( inc_x1 != 1 ||
	     inc_x2 != 1 ||
	     inc_y  != 1 ) bl1_abort();

	n_pre = 0;
	if ( ( unsigned long ) y % 16 != 0 )
	{
		if ( ( unsigned long ) x1 % 16 == 0 ||
		     ( unsigned long ) x2 % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	chi1 = x1;
	chi2 = x2;
	psi1 = y;

	if ( n_pre == 1 )
	{
		double   alpha1_c = *alpha1;
		double   alpha2_c = *alpha2;
		double   chi11_c = *chi1;
		double   chi12_c = *chi2;
		double   temp1;

		// psi1 = psi1 + alpha1 * chi11 + alpha2 * chi12;
		temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;
		*psi1 = *psi1 + temp1;

		chi1 += inc_x1;
		chi2 += inc_x2;
		psi1 += inc_y;
	}

	a1v.v = _mm_loaddup_pd( ( double* )alpha1 );
	a2v.v = _mm_loaddup_pd( ( double* )alpha2 );

	for ( i = 0; i < n_run; ++i )
	{
		x11v.v = _mm_load_pd( ( double* )chi1 );
		x12v.v = _mm_load_pd( ( double* )chi2 );
		y1v.v  = _mm_load_pd( ( double* )psi1 );

		x21v.v = _mm_load_pd( ( double* )(chi1 + 2) );
		x22v.v = _mm_load_pd( ( double* )(chi2 + 2) );
		y2v.v  = _mm_load_pd( ( double* )(psi1 + 2) );

		y1v.v += a1v.v * x11v.v + a2v.v * x12v.v;
		y2v.v += a1v.v * x21v.v + a2v.v * x22v.v;

		_mm_store_pd( ( double* )psi1, y1v.v );
		_mm_store_pd( ( double* )(psi1 + 2), y2v.v );

		//chi1 += step_x1;
		//chi2 += step_x2;
		//psi1 += step_y;
		chi1 += 4;
		chi2 += 4;
		psi1 += 4;
	}

	if ( n_left > 0 )
	{
		double   alpha1_c = *alpha1;
		double   alpha2_c = *alpha2;

		for ( i = 0; i < n_left; ++i )
		{
			double   chi11_c = *chi1;
			double   chi12_c = *chi2;
			double   psi1_c  = *psi1;
			double   temp1;

			temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;
			*psi1 = psi1_c + temp1;

			chi1 += inc_x1;
			chi2 += inc_x2;
			psi1 += inc_y;
		}
	}
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict chi2;
	double*   restrict psi1;
	double    alpha1_c;
	double    alpha2_c;
	double    temp1;
	double    temp2;
	int       i;

	int       n_run  = n / 2;
	int       n_left = n % 2;
	int       twoinc_x1 = 2*inc_x1;
	int       twoinc_x2 = 2*inc_x2;
	int       twoinc_y  = 2*inc_y;

	chi1 = x1;
	chi2 = x2;
	psi1 = y;

	alpha1_c = *alpha1;
	alpha2_c = *alpha2;

	for ( i = 0; i < n_run; ++i )
	{
		double   chi11_c = *chi1;
		double   chi21_c = *(chi1 + inc_x1);
		double   chi12_c = *chi2;
		double   chi22_c = *(chi2 + inc_x2);
		double   psi1_c  = *psi1;
		double   psi2_c  = *(psi1 + inc_y);

		// psi1 = psi1 + alpha1 * chi11 + alpha2 * chi12;
		// psi2 = psi2 + alpha1 * chi21 + alpha2 * chi22;
		temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;
		temp2 = alpha1_c * chi21_c + alpha2_c * chi22_c;

		*psi1           = psi1_c + temp1;
		*(psi1 + inc_y) = psi2_c + temp2;

		chi1 += twoinc_x1;
		chi2 += twoinc_x2;
		psi1 += twoinc_y;
	}

	if ( n_left == 1 )
	{
		double   chi11_c = *chi1;
		double   chi12_c = *chi2;

		// psi1 = psi1 + alpha1 * chi11 + alpha2 * chi12;
		temp1 = alpha1_c * chi11_c + alpha2_c * chi12_c;

		*psi1 = *psi1 + temp1;
	}
}
#endif


void bl1_caxpyv2b( int       n,
                   scomplex* alpha1,
                   scomplex* alpha2,
                   scomplex* x1, int inc_x1,
                   scomplex* x2, int inc_x2,
                   scomplex* y,  int inc_y )
{
	bl1_abort();
}


void bl1_zaxpyv2b( int       n,
                   dcomplex* alpha1,
                   dcomplex* alpha2,
                   dcomplex* x1, int inc_x1,
                   dcomplex* x2, int inc_x2,
                   dcomplex* y,  int inc_y )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict chi2;
	dcomplex* restrict psi1;
	int       i;
	v2df_t    alpha1v, alpha1rv;
	v2df_t    alpha2v, alpha2rv;
	v2df_t    x11v, x12v;
	v2df_t    t1v, y1v;
	v2df_t    acbc, bdad;

	chi1 = x1;
	chi2 = x2;
	psi1 = y;

	alpha1v.v  = _mm_load_pd( ( double* )alpha1 );
	alpha2v.v  = _mm_load_pd( ( double* )alpha2 );
	alpha1rv.v = _mm_shuffle_pd( alpha1v.v, alpha1v.v, _MM_SHUFFLE2 (0,1) );
	alpha2rv.v = _mm_shuffle_pd( alpha2v.v, alpha2v.v, _MM_SHUFFLE2 (0,1) );

	if ( inc_x1 == 1 &&
	     inc_x2 == 1 &&
	     inc_y  == 1 )
	{
		for ( i = 0; i < n; ++i )
		{
			x11v.v = _mm_load_pd( ( double* )chi1 );
			x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
			x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );
			acbc.v = alpha1v.v  * x11v.v;
			bdad.v = alpha1rv.v * x12v.v;
			t1v.v = _mm_addsub_pd( acbc.v, bdad.v );

			x11v.v = _mm_load_pd( ( double* )chi2 );
			x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
			x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );
			acbc.v = alpha2v.v  * x11v.v;
			bdad.v = alpha2rv.v * x12v.v;
			t1v.v = t1v.v + _mm_addsub_pd( acbc.v, bdad.v );

			y1v.v = _mm_load_pd( ( double* )psi1 );
			y1v.v = y1v.v + t1v.v;
			_mm_store_pd( ( double* )psi1, y1v.v );

			chi1 += 1;
			chi2 += 1;
			psi1 += 1;
		}
	}
	else
	{
		for ( i = 0; i < n; ++i )
		{
			x11v.v = _mm_load_pd( ( double* )chi1 );
			x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
			x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );
			acbc.v = alpha1v.v  * x11v.v;
			bdad.v = alpha1rv.v * x12v.v;
			t1v.v = _mm_addsub_pd( acbc.v, bdad.v );

			x11v.v = _mm_load_pd( ( double* )chi2 );
			x12v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (1,1) );
			x11v.v = _mm_shuffle_pd( x11v.v, x11v.v, _MM_SHUFFLE2 (0,0) );
			acbc.v = alpha2v.v  * x11v.v;
			bdad.v = alpha2rv.v * x12v.v;
			t1v.v = t1v.v + _mm_addsub_pd( acbc.v, bdad.v );

			y1v.v = _mm_load_pd( ( double* )psi1 );
			y1v.v = y1v.v + t1v.v;
			_mm_store_pd( ( double* )psi1, y1v.v );

			chi1 += inc_x1;
			chi2 += inc_x2;
			psi1 += inc_y;
		}
	}
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict chi2;
	dcomplex* restrict psi1;
	dcomplex  alpha1_c;
	dcomplex  alpha2_c;
	dcomplex  temp;
	int       i;

	chi1 = x1;
	chi2 = x2;
	psi1 = y;

	alpha1_c = *alpha1;
	alpha2_c = *alpha2;

	for ( i = 0; i < n; ++i )
	{
		dcomplex chi1_c = *chi1;
		dcomplex chi2_c = *chi2;

		temp.real = 0.0;
		temp.imag = 0.0;

		// psi1 = psi1 + alpha1 * chi1;
		temp.real += alpha1_c.real * chi1_c.real - alpha1_c.imag * chi1_c.imag;
		temp.imag += alpha1_c.real * chi1_c.imag + alpha1_c.imag * chi1_c.real;

		// psi1 = psi1 + alpha2 * chi2;
		temp.real += alpha2_c.real * chi2_c.real - alpha2_c.imag * chi2_c.imag;
		temp.imag += alpha2_c.real * chi2_c.imag + alpha2_c.imag * chi2_c.real;

		psi1->real = psi1->real + temp.real;
		psi1->imag = psi1->imag + temp.imag;

		chi1 += inc_x1;
		chi2 += inc_x2;
		psi1 += inc_y;
	}
}
#endif

