
#include "blis1.h"

/*
   Effective computation:

     a   = a + beta * u + gamma * z;
     rho = conj(a) * x;
     w   = w + kappa * a;
*/

void bl1_saxpyv2bdotaxpy( int       n,
                          float*    beta,
                          float*    u, int inc_u,
                          float*    gamma,
                          float*    z, int inc_z,
                          float*    a, int inc_a,
                          float*    x, int inc_x,
                          float*    kappa,
                          float*    rho,
                          float*    w, int inc_w )
{
	bl1_abort();
}


void bl1_daxpyv2bdotaxpy( int       n,
                          double*   beta,
                          double*   u, int inc_u,
                          double*   gamma,
                          double*   z, int inc_z,
                          double*   a, int inc_a,
                          double*   x, int inc_x,
                          double*   kappa,
                          double*   rho,
                          double*   w, int inc_w )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	double*   restrict upsilon1;
	double*   restrict zeta1;
	double*   restrict alpha1;
	double*   restrict chi1;
	double*   restrict omega1;
	double             rho_c;
	int                i;
	v2df_t             b1v, g1v, k1v;
	v2df_t             rhov;
	v2df_t             u1v, z1v, a1v;
	v2df_t             u2v, z2v, a2v;
	v2df_t             x1v, w1v;
	v2df_t             x2v, w2v;

	int       n_pre;
	int       n_run;
	int       n_left;

	n_pre = 0;
	if ( ( unsigned long ) a % 16 != 0 )
	{
		if ( ( unsigned long ) u % 16 == 0 ||
		     ( unsigned long ) z % 16 == 0 ||
		     ( unsigned long ) x % 16 == 0 ||
		     ( unsigned long ) w % 16 == 0 ) bl1_abort();

		n_pre = 1;
	}

	n_run       = ( n - n_pre ) / 4;
	n_left      = ( n - n_pre ) % 4;

	upsilon1 = u;
	zeta1    = z;
	alpha1   = a;
	chi1     = x;
	omega1   = w;


	rho_c   = 0.0;

	if ( n_pre == 1 )
	{
		double   beta_c     = *beta;
		double   gamma_c    = *gamma;
		double   kappa_c    = *kappa;

		double   upsilon1_c = *upsilon1;
		double   zeta1_c    = *zeta1;
		double   alpha1_c   = *alpha1;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		alpha1_c += beta_c * upsilon1_c + gamma_c * zeta1_c;
		rho_c += alpha1_c * chi1_c;
		omega1_c += kappa_c * alpha1_c;

		*alpha1 = alpha1_c;
		*omega1 = omega1_c;

		upsilon1 += inc_u;
		zeta1    += inc_z;
		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	b1v.v = _mm_loaddup_pd( ( double* )beta );
	g1v.v = _mm_loaddup_pd( ( double* )gamma );
	k1v.v = _mm_loaddup_pd( ( double* )kappa );

	rhov.v = _mm_setzero_pd();

	for ( i = 0; i < n_run; ++i )
	{
		u1v.v = _mm_load_pd( ( double* )upsilon1 );
		z1v.v = _mm_load_pd( ( double* )zeta1 );
		a1v.v = _mm_load_pd( ( double* )alpha1 );

		a1v.v += b1v.v * u1v.v + g1v.v * z1v.v;

		u2v.v = _mm_load_pd( ( double* )(upsilon1 + 2) );
		z2v.v = _mm_load_pd( ( double* )(zeta1 + 2) );
		a2v.v = _mm_load_pd( ( double* )(alpha1 + 2) );

		a2v.v += b1v.v * u2v.v + g1v.v * z2v.v;

		x1v.v = _mm_load_pd( ( double* )chi1 );
		x2v.v = _mm_load_pd( ( double* )(chi1 + 2) );

		w1v.v = _mm_load_pd( ( double* )omega1 );
		w2v.v = _mm_load_pd( ( double* )(omega1 + 2) );

		rhov.v += a1v.v * x1v.v;
		rhov.v += a2v.v * x2v.v;

		w1v.v += k1v.v * a1v.v;
		w2v.v += k1v.v * a2v.v;

		_mm_store_pd( ( double* )alpha1, a1v.v );
		_mm_store_pd( ( double* )(alpha1 + 2), a2v.v );

		_mm_store_pd( ( double* )omega1, w1v.v );
		_mm_store_pd( ( double* )(omega1 + 2), w2v.v );


		upsilon1 += 4;
		zeta1    += 4;
		alpha1   += 4;
		chi1     += 4;
		omega1   += 4;
	}

	rho_c += rhov.d[0] + rhov.d[1];

	if ( n_left > 0 )
	{
		double beta_c  = *beta;
		double gamma_c = *gamma;
		double kappa_c = *kappa;

		for ( i = 0; i < n_left; ++i )
		{
			double   upsilon1_c = *upsilon1;
			double   zeta1_c    = *zeta1;
			double   alpha1_c   = *alpha1;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			alpha1_c += beta_c * upsilon1_c + gamma_c * zeta1_c;
			rho_c += alpha1_c * chi1_c;
			omega1_c += kappa_c * alpha1_c;

			*alpha1 = alpha1_c;
			*omega1 = omega1_c;

			upsilon1 += inc_u;
			zeta1    += inc_z;
			alpha1   += inc_a;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	*rho = rho_c;
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	double*   restrict upsilon1;
	double*   restrict zeta1;
	double*   restrict alpha1;
	double*   restrict chi1;
	double*   restrict omega1;
	double             beta_c;
	double             gamma_c;
	double             kappa_c;
	double             rho_c;
	int                i;

	int       n_pre;
	int       n_run;
	int       n_left;

	n_pre = 0;
	//if ( ( unsigned long ) a % 16 != 0 )
	//{
	//	if ( ( unsigned long ) u % 16 == 0 ||
	//	     ( unsigned long ) z % 16 == 0 ||
	//	     ( unsigned long ) x % 16 == 0 ||
	//	     ( unsigned long ) w % 16 == 0 ) bl1_abort();
	//
	//	n_pre = 1;
	//}

	n_run       = ( n - n_pre ) / 2;
	n_left      = ( n - n_pre ) % 2;

	upsilon1 = u;
	zeta1    = z;
	alpha1   = a;
	chi1     = x;
	omega1   = w;

	beta_c  = *beta;
	gamma_c = *gamma;
	kappa_c = *kappa;

	rho_c = 0.0;

	if ( n_pre == 1 )
	{
		double   upsilon1_c = *upsilon1;
		double   zeta1_c    = *zeta1;
		double   alpha1_c   = *alpha1;
		double   chi1_c     = *chi1;
		double   omega1_c   = *omega1;

		alpha1_c += beta_c * upsilon1_c + gamma_c * zeta1_c;
		rho_c += alpha1_c * chi1_c;
		omega1_c += kappa_c * alpha1_c;

		*alpha1 = alpha1_c;
		*omega1 = omega1_c;

		upsilon1 += inc_u;
		zeta1    += inc_z;
		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	for ( i = 0; i < n_run; ++i )
	{
		double   upsilon1_c = *upsilon1;
		double   upsilon2_c = *(upsilon1 + 1);
		double   zeta1_c    = *zeta1;
		double   zeta2_c    = *(zeta1 + 1);
		double   alpha1_c   = *alpha1;
		double   alpha2_c   = *(alpha1 + 1);
		double   chi1_c     = *chi1;
		double   chi2_c     = *(chi1 + 1);
		double   omega1_c   = *omega1;
		double   omega2_c   = *(omega1 + 1);

		// alpha1 += beta * upsilon1 + gamma * zeta1;
		alpha1_c += beta_c * upsilon1_c + gamma_c * zeta1_c;
		alpha2_c += beta_c * upsilon2_c + gamma_c * zeta2_c;

		// rho += conj(alpha1) * chi1 +
		//        conj(alpha2) * chi2;
		rho_c += alpha1_c * chi1_c + alpha2_c * chi2_c;

		// omega1 += kappa * alpha1;
		omega1_c += kappa_c * alpha1_c;
		omega2_c += kappa_c * alpha2_c;

		*alpha1       = alpha1_c;
		*(alpha1 + 1) = alpha2_c;
		*omega1       = omega1_c;
		*(omega1 + 1) = omega2_c;

		upsilon1 += 2*inc_u;
		zeta1    += 2*inc_z;
		alpha1   += 2*inc_a;
		chi1     += 2*inc_x;
		omega1   += 2*inc_w;
	}

	if ( n_left > 0 )
	{

		for ( i = 0; i < n_left; ++i )
		{
			double   upsilon1_c = *upsilon1;
			double   zeta1_c    = *zeta1;
			double   alpha1_c   = *alpha1;
			double   chi1_c     = *chi1;
			double   omega1_c   = *omega1;

			alpha1_c += beta_c * upsilon1_c + gamma_c * zeta1_c;
			rho_c += alpha1_c * chi1_c;
			omega1_c += kappa_c * alpha1_c;

			*alpha1 = alpha1_c;
			*omega1 = omega1_c;

			upsilon1 += inc_u;
			zeta1    += inc_z;
			alpha1   += inc_a;
			chi1     += inc_x;
			omega1   += inc_w;
		}
	}

	*rho = rho_c;
}
#endif


void bl1_caxpyv2bdotaxpy( int       n,
                          scomplex* beta,
                          scomplex* u, int inc_u,
                          scomplex* gamma,
                          scomplex* z, int inc_z,
                          scomplex* a, int inc_a,
                          scomplex* x, int inc_x,
                          scomplex* kappa,
                          scomplex* rho,
                          scomplex* w, int inc_w )
{
	bl1_abort();
}


void bl1_zaxpyv2bdotaxpy( int       n,
                          dcomplex* beta,
                          dcomplex* u, int inc_u,
                          dcomplex* gamma,
                          dcomplex* z, int inc_z,
                          dcomplex* a, int inc_a,
                          dcomplex* x, int inc_x,
                          dcomplex* kappa,
                          dcomplex* rho,
                          dcomplex* w, int inc_w )
#if BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_SSE_INTRINSICS
{
	dcomplex* restrict upsilon1;
	dcomplex* restrict zeta1;
	dcomplex* restrict alpha1;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	int                i;

	//v2df_t    beta1v, beta1rv;
	//v2df_t    gamma1v, gamma1rv;
	//v2df_t    kappa1v, kappa1rv;
	v2df_t    rho1v;
	//v2df_t    u11v, u12v;
	//v2df_t    z11v, z12v;
	v2df_t    a11v, a12v;
	v2df_t    x1v, x1rv;
	v2df_t    w1v;
	v2df_t    acbc, bdad;
	v2df_t    adac, bcbd;

	v2df_t    a1v, a1rv;
	v2df_t    u1v, u1rv;
	v2df_t    z1v, z1rv;
	v2df_t    beta11v, gamma11v, kappa11v;
	v2df_t    beta12v, gamma12v, kappa12v;

	upsilon1 = u;
	zeta1    = z;
	alpha1   = a;
	chi1     = x;
	omega1   = w;

	if ( inc_u != 1 || 
	     inc_z != 1 ||
	     inc_a != 1 ||
	     inc_x != 1 ||
	     inc_w != 1 ) bl1_abort();


	beta11v.v  = _mm_loaddup_pd( ( double* )&(beta->real) );
	beta12v.v  = _mm_loaddup_pd( ( double* )&(beta->imag) );
	gamma11v.v = _mm_loaddup_pd( ( double* )&(gamma->real) );
	gamma12v.v = _mm_loaddup_pd( ( double* )&(gamma->imag) );
	kappa11v.v = _mm_loaddup_pd( ( double* )&(kappa->real) );
	kappa12v.v = _mm_loaddup_pd( ( double* )&(kappa->imag) );

	rho1v.v = _mm_setzero_pd();

	for ( i = 0; i < n; ++i )
	{
		//alpha_c = *alpha1;
		a1v.v  = _mm_load_pd( ( double* )alpha1 );

		//alpha1_c.real += beta_c.real * upsilon1_c.real - beta_c.imag * upsilon1_c.imag;
		//alpha1_c.imag += beta_c.real * upsilon1_c.imag + beta_c.imag * upsilon1_c.real;
		u1v.v  = _mm_load_pd( ( double* )upsilon1 );
		u1rv.v = _mm_shuffle_pd( u1v.v, u1v.v, _MM_SHUFFLE2 (0,1) );
		acbc.v = beta11v.v * u1v.v;
		bdad.v = beta12v.v * u1rv.v;
		a1v.v += _mm_addsub_pd( acbc.v, bdad.v );

		//alpha1_c.real += gamma_c.real * zeta1_c.real - gamma_c.imag * zeta1_c.imag;
		//alpha1_c.imag += gamma_c.real * zeta1_c.imag + gamma_c.imag * zeta1_c.real;
		z1v.v  = _mm_load_pd( ( double* )zeta1 );
		z1rv.v = _mm_shuffle_pd( z1v.v, z1v.v, _MM_SHUFFLE2 (0,1) );
		acbc.v = gamma11v.v * z1v.v;
		bdad.v = gamma12v.v * z1rv.v;
		a1v.v += _mm_addsub_pd( acbc.v, bdad.v );

		//*alpha1 = alpha1_c;
		_mm_store_pd( ( double* )alpha1, a1v.v );

		//rho_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		//rho_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;
		x1v.v  = _mm_load_pd( ( double* )chi1 );
		x1rv.v = _mm_shuffle_pd( x1v.v, x1v.v, _MM_SHUFFLE2 (0,1) );
		a11v.v = a1v.v;
		a12v.v = _mm_shuffle_pd( a11v.v, a11v.v, _MM_SHUFFLE2 (1,1) );
		a11v.v = _mm_shuffle_pd( a11v.v, a11v.v, _MM_SHUFFLE2 (0,0) );
		adac.v = a11v.v * x1rv.v;
		bcbd.v = a12v.v * x1v.v;
		rho1v.v = rho1v.v + _mm_addsub_pd( adac.v, bcbd.v );

		//omega_c = *omega1;
		w1v.v  = _mm_load_pd( ( double* )omega1 );

		//omega1_c.real += kappa_c.real * alpha1_c.real - kappa_c.imag * alpha1_c.imag;
		//omega1_c.imag += kappa_c.real * alpha1_c.imag + kappa_c.imag * alpha1_c.real;
		a1rv.v = _mm_shuffle_pd( a1v.v, a1v.v, _MM_SHUFFLE2 (0,1) );
		acbc.v = kappa11v.v * a1v.v;
		bdad.v = kappa12v.v * a1rv.v;
		w1v.v += _mm_addsub_pd( acbc.v, bdad.v );

		// *omega1 = omega1_c;
		_mm_store_pd( ( double* )omega1, w1v.v );


		upsilon1 += 1;
		zeta1    += 1;
		alpha1   += 1;
		chi1     += 1;
		omega1   += 1;
	}

	rho1v.v = _mm_shuffle_pd( rho1v.v, rho1v.v, _MM_SHUFFLE2 (0,1) );

	//rho->real = rho_c.real;
	//rho->imag = rho_c.imag;
	_mm_store_pd( ( double* )rho, rho1v.v );
}
#elif BLIS1_VECTOR_INTRINSIC_TYPE == BLIS1_NO_INTRINSICS
{
	dcomplex* restrict upsilon1;
	dcomplex* restrict zeta1;
	dcomplex* restrict alpha1;
	dcomplex* restrict chi1;
	dcomplex* restrict omega1;
	dcomplex           beta_c;
	dcomplex           gamma_c;
	dcomplex           kappa_c;
	dcomplex           rho_c;
	int                i;

	upsilon1 = u;
	zeta1    = z;
	alpha1   = a;
	chi1     = x;
	omega1   = w;

	rho_c.real = 0.0;
	rho_c.imag = 0.0;

	beta_c  = *beta;
	gamma_c = *gamma;
	kappa_c = *kappa;

	for ( i = 0; i < n; ++i )
	{
		dcomplex upsilon1_c = *upsilon1;
		dcomplex zeta1_c    = *zeta1;
		dcomplex alpha1_c   = *alpha1;
		dcomplex chi1_c     = *chi1;
		dcomplex omega1_c   = *omega1;

		// alpha1 += beta * upsilon1;
		alpha1_c.real += beta_c.real * upsilon1_c.real - beta_c.imag * upsilon1_c.imag;
		alpha1_c.imag += beta_c.real * upsilon1_c.imag + beta_c.imag * upsilon1_c.real;
		
		// alpha1 += gamma * zeta1;
		alpha1_c.real += gamma_c.real * zeta1_c.real - gamma_c.imag * zeta1_c.imag;
		alpha1_c.imag += gamma_c.real * zeta1_c.imag + gamma_c.imag * zeta1_c.real;
		
		// rho += conj(alpha1) * chi1;
		rho_c.real += alpha1_c.real * chi1_c.real - -alpha1_c.imag * chi1_c.imag;
		rho_c.imag += alpha1_c.real * chi1_c.imag + -alpha1_c.imag * chi1_c.real;

		// omega1 += kappa * alpha1;
		omega1_c.real += kappa_c.real * alpha1_c.real - kappa_c.imag * alpha1_c.imag;
		omega1_c.imag += kappa_c.real * alpha1_c.imag + kappa_c.imag * alpha1_c.real;

		*alpha1 = alpha1_c;
		*omega1 = omega1_c;

		upsilon1 += inc_u;
		zeta1    += inc_z;
		alpha1   += inc_a;
		chi1     += inc_x;
		omega1   += inc_w;
	}

	rho->real = rho_c.real;
	rho->imag = rho_c.imag;
}
#endif

