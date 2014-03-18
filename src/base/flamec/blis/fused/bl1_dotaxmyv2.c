
#include "FLAME.h"

/*
   Effective computation:

     rho = conj(x) * u;
     y   = y - alpha * x;
     z   = z - beta  * x;
*/

void bl1_sdotaxmyv2( int       n,
                     float*    alpha,
                     float*    beta,
                     float*    x, int inc_x,
                     float*    u, int inc_u,
                     float*    rho,
                     float*    y, int inc_y,
                     float*    z, int inc_z )
{
	bl1_abort();
}


void bl1_ddotaxmyv2( int       n,
                     double*   alpha,
                     double*   beta,
                     double*   x, int inc_x,
                     double*   u, int inc_u,
                     double*   rho,
                     double*   y, int inc_y,
                     double*   z, int inc_z )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict upsilon1;
	double*   restrict psi1;
	double*   restrict zeta1;
	double    rho_c;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	v2df_t    a1v, b1v;
	v2df_t    rho1v;
	v2df_t    x1v, u1v, y1v, z1v;

	if ( inc_x != 1 ||
	     inc_u != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	if ( ( unsigned long ) z % 16 != 0 )
	{
		if ( ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) u % 16 == 0 ||
		     ( unsigned long ) y % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	chi1     = x;
	upsilon1 = u;
	psi1     = y;
	zeta1    = z;

	rho_c   = 0.0;

	if ( n_pre == 1 )
	{
		double   alpha_c = *alpha;
		double   beta_c  = *beta;
		double   chi1_c    = *chi1;
		double   upsilon_c = *upsilon1;

		rho_c  += chi1_c * upsilon_c;
		*psi1  -= alpha_c * chi1_c;
		*zeta1 -= beta_c  * chi1_c;

		chi1     += inc_x;
		upsilon1 += inc_u;
		psi1     += inc_y;
		zeta1    += inc_z;
	}

	a1v.v = _mm_loaddup_pd( ( double* )alpha );
	b1v.v = _mm_loaddup_pd( ( double* )beta );

	rho1v.v = _mm_setzero_pd();

	for ( i = 0; i < n_run; ++i )
	{
		x1v.v = _mm_load_pd( ( double* )chi1 );
		u1v.v = _mm_load_pd( ( double* )upsilon1 );
		y1v.v = _mm_load_pd( ( double* )psi1 );
		z1v.v = _mm_load_pd( ( double* )zeta1 );

		rho1v.v += x1v.v * u1v.v;
		y1v.v   -= a1v.v * x1v.v;
		z1v.v   -= b1v.v * x1v.v;

		_mm_store_pd( ( double* )psi1,  y1v.v );
		_mm_store_pd( ( double* )zeta1, z1v.v );

		chi1     += 2;
		upsilon1 += 2;
		psi1     += 2;
		zeta1    += 2;
	}

	rho_c += rho1v.d[0] + rho1v.d[1];

	if ( n_left > 0 )
	{
		double   alpha_c = *alpha;
		double   beta_c  = *beta;

		for( i = 0; i < n_left; ++i )
		{
			double   chi1_c    = *chi1;
			double   upsilon_c = *upsilon1;

			rho_c  += chi1_c * upsilon_c;
			*psi1  -= alpha_c * chi1_c;
			*zeta1 -= beta_c  * chi1_c;

			chi1     += inc_x;
			upsilon1 += inc_u;
			psi1     += inc_y;
			zeta1    += inc_z;
		}
	}

	*rho = rho_c;
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict chi1;
	double*   restrict upsilon1;
	double*   restrict psi1;
	double*   restrict zeta1;
	double    alpha_c;
	double    beta_c;
	double    rho_c;
	int       i;

	int       n_pre;
	int       n_run;
	int       n_left;

	if ( inc_x != 1 ||
	     inc_u != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	n_pre = 0;
	//if ( ( unsigned long ) z % 16 != 0 )
	//{
	//	if ( ( unsigned long ) x % 16 == 0 ||
	//	     ( unsigned long ) u % 16 == 0 ||
	//	     ( unsigned long ) y % 16 == 0 ) bl1_abort();
	//
	//	n_pre = 1;
	//}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	chi1     = x;
	upsilon1 = u;
	psi1     = y;
	zeta1    = z;

	alpha_c = *alpha;
	beta_c  = *beta;

	rho_c   = 0.0;

	if ( n_pre == 1 )
	{
		double   chi1_c    = *chi1;
		double   upsilon_c = *upsilon1;

		rho_c  += chi1_c * upsilon_c;
		*psi1  -= alpha_c * chi1_c;
		*zeta1 -= beta_c  * chi1_c;

		chi1     += inc_x;
		upsilon1 += inc_u;
		psi1     += inc_y;
		zeta1    += inc_z;
	}

	for ( i = 0; i < n_run; ++i )
	{
		double   chi1_c     = *chi1;
		double   chi2_c     = *(chi1 + 1);
		double   upsilon1_c = *upsilon1;
		double   upsilon2_c = *(upsilon1 + 1);
		double   psi1_c     = *psi1;
		double   psi2_c     = *(psi1 + 1);
		double   zeta1_c    = *zeta1;
		double   zeta2_c    = *(zeta1 + 1);

		rho_c   += chi1_c * upsilon1_c;
		rho_c   += chi2_c * upsilon2_c;

		psi1_c  -= alpha_c * chi1_c;
		psi2_c  -= alpha_c * chi2_c;

		zeta1_c -= beta_c  * chi1_c;
		zeta2_c -= beta_c  * chi2_c;

		*psi1        = psi1_c;
		*(psi1 + 1)  = psi2_c;
		*zeta1       = zeta1_c;
		*(zeta1 + 1) = zeta2_c;

		chi1     += 2*inc_x;
		upsilon1 += 2*inc_u;
		psi1     += 2*inc_y;
		zeta1    += 2*inc_z;
	}

	if ( n_left > 0 )
	{
		for( i = 0; i < n_left; ++i )
		{
			double   chi1_c    = *chi1;
			double   upsilon_c = *upsilon1;

			rho_c  += chi1_c * upsilon_c;
			*psi1  -= alpha_c * chi1_c;
			*zeta1 -= beta_c  * chi1_c;

			chi1     += inc_x;
			upsilon1 += inc_u;
			psi1     += inc_y;
			zeta1    += inc_z;
		}
	}

	*rho = rho_c;
}
#endif


void bl1_cdotaxmyv2( int       n,
                     scomplex* alpha,
                     scomplex* beta,
                     scomplex* x, int inc_x,
                     scomplex* u, int inc_u,
                     scomplex* rho,
                     scomplex* y, int inc_y,
                     scomplex* z, int inc_z )
{
	bl1_abort();
}


void bl1_zdotaxmyv2( int       n,
                     dcomplex* alpha,
                     dcomplex* beta,
                     dcomplex* x, int inc_x,
                     dcomplex* u, int inc_u,
                     dcomplex* rho,
                     dcomplex* y, int inc_y,
                     dcomplex* z, int inc_z )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict upsilon1;
	dcomplex* restrict psi1;
	dcomplex* restrict zeta1;
	int       i;

	v2df_t alpha11v, alpha12v;
	v2df_t beta11v,  beta12v;
	v2df_t rho1v;
	v2df_t x1v, x1rv;
	v2df_t y1v;
	v2df_t z1v;
	v2df_t u11v, u12v;
	v2df_t acad, bdbc;
	v2df_t bcac, adbd;

	if ( inc_x != 1 ||
	     inc_u != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	chi1     = x;
	upsilon1 = u;
	psi1     = y;
	zeta1    = z;

	//rho_c.real = 0.0;
	//rho_c.imag = 0.0;
	rho1v.v = _mm_setzero_pd();

	//alpha_c = *alpha;
	//beta_c  = *beta;
	alpha11v.v  = _mm_loaddup_pd( ( double* )&(alpha->real) );
	alpha12v.v  = _mm_loaddup_pd( ( double* )&(alpha->imag) );
	beta11v.v   = _mm_loaddup_pd( ( double* )&(beta->real) );
	beta12v.v   = _mm_loaddup_pd( ( double* )&(beta->imag) );

	for ( i = 0; i < n; ++i )
	{
		//dcomplex chi1_c     = *chi1;
		x1v.v = _mm_load_pd( ( double* )chi1 );

		//psi1->real  -= alpha_c.real * chi1_c.real - alpha_c.imag * chi1_c.imag;
		//psi1->imag  -= alpha_c.real * chi1_c.imag + alpha_c.imag * chi1_c.real;
		x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
		acad.v = alpha11v.v * x1v.v;
		bdbc.v = alpha12v.v * x1rv.v;
		y1v.v = _mm_load_pd( ( double* )psi1 );
		y1v.v = y1v.v - _mm_addsub_pd( acad.v, bdbc.v );
		_mm_store_pd( ( double* )psi1, y1v.v );

		//zeta1->real -= beta_c.real  * chi1_c.real - beta_c.imag  * chi1_c.imag;
		//zeta1->imag -= beta_c.real  * chi1_c.imag + beta_c.imag  * chi1_c.real;
		x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
		acad.v = beta11v.v * x1v.v;
		bdbc.v = beta12v.v * x1rv.v;
		z1v.v = _mm_load_pd( ( double* )zeta1 );
		z1v.v = z1v.v - _mm_addsub_pd( acad.v, bdbc.v );
		_mm_store_pd( ( double* )zeta1, z1v.v );

		//rho_c.real = chi1_c.real * upsilon1_c.real - -chi1_c.imag * upsilon1_c.imag;
		//rho_c.imag = chi1_c.real * upsilon1_c.imag + -chi1_c.imag * upsilon1_c.real;
		x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
		u11v.v = _mm_loaddup_pd( ( double* )&(upsilon1->real) );
		u12v.v = _mm_loaddup_pd( ( double* )&(upsilon1->imag) );
		bcac.v = x1rv.v * u11v.v;
		adbd.v = x1v.v  * u12v.v;
		rho1v.v = rho1v.v + _mm_addsub_pd( bcac.v, adbd.v );

		chi1     += 1;
		upsilon1 += 1;
		psi1     += 1;
		zeta1    += 1;
	}

	rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );

	rho1v.d[1] = -rho1v.d[1];

	_mm_store_pd( ( double* )rho, rho1v.v );
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict chi1;
	dcomplex* restrict upsilon1;
	dcomplex* restrict psi1;
	dcomplex* restrict zeta1;
	dcomplex  alpha_c;
	dcomplex  beta_c;
	dcomplex  rho_c;
	int       i;

	if ( inc_x != 1 ||
	     inc_u != 1 ||
	     inc_y != 1 ||
	     inc_z != 1 ) bl1_abort();

	chi1     = x;
	upsilon1 = u;
	psi1     = y;
	zeta1    = z;

	rho_c.real = 0.0;
	rho_c.imag = 0.0;

	alpha_c = *alpha;
	beta_c  = *beta;

	for ( i = 0; i < n; ++i )
	{
		dcomplex chi1_c     = *chi1;
		dcomplex upsilon1_c = *upsilon1;

		// rho += conj(chi1) * upsilon1;
		rho_c.real += chi1_c.real * upsilon1_c.real - -chi1_c.imag * upsilon1_c.imag;
		rho_c.imag += chi1_c.real * upsilon1_c.imag + -chi1_c.imag * upsilon1_c.real;

		// psi1  = psi1  - alpha * chi1;
		psi1->real  -= alpha_c.real * chi1_c.real - alpha_c.imag * chi1_c.imag;
		psi1->imag  -= alpha_c.real * chi1_c.imag + alpha_c.imag * chi1_c.real;

		// zeta1 = zeta1 - beta  * chi1;
		zeta1->real -= beta_c.real  * chi1_c.real - beta_c.imag  * chi1_c.imag;
		zeta1->imag -= beta_c.real  * chi1_c.imag + beta_c.imag  * chi1_c.real;

		chi1     += inc_x;
		upsilon1 += inc_u;
		psi1     += inc_y;
		zeta1    += inc_z;
	}
	
	*rho = rho_c;
}
#endif

