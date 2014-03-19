/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#if FLA_VECTOR_INTRINSIC_TYPE == FLA_NO_INTRINSICS

#define MAC_Apply_G_mx3_ass MAC_Apply_G_mx3_ops
#define MAC_Apply_G_mx3_asd MAC_Apply_G_mx3_opd
#define MAC_Apply_G_mx3_asc MAC_Apply_G_mx3_opc
#define MAC_Apply_G_mx3_asz MAC_Apply_G_mx3_opz

#elif FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS

#define MAC_Apply_G_mx3_ass( m_A, \
                             gamma12, \
                             sigma12, \
                             gamma23, \
                             sigma23, \
                             a1, inc_a1, \
                             a2, inc_a2, \
                             a3, inc_a3 ) \
{\
	int              n_iter32 = m_A / ( 4 * 8 ); \
	int              n_left32 = m_A % ( 4 * 8 ); \
	int              n_iter4  = n_left32 / ( 4 * 1 ); \
	int              n_left   = n_left32 % ( 4 * 1 ); \
	int              i; \
\
	const int        step_a1 = inc_a1 * 4; \
	const int        step_a2 = inc_a1 * 4; \
	const int        step_a3 = inc_a1 * 4; \
\
	float* restrict alpha1 = a1; \
	float* restrict alpha2 = a2; \
	float* restrict alpha3 = a3; \
\
	v4sf_t    a1v, a2v, a3v; \
	v4sf_t    g12v, s12v; \
	v4sf_t    g23v, s23v; \
	v4sf_t    t1v, t2v; \
\
	g12v.v = _mm_load1_ps( gamma12 ); \
	s12v.v = _mm_load1_ps( sigma12 ); \
	g23v.v = _mm_load1_ps( gamma23 ); \
	s23v.v = _mm_load1_ps( sigma23 ); \
\
	for ( i = 0; i < n_iter32; ++i ) \
	{ \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_iter4; ++i ) \
	{ \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
		float ga12 = *gamma12; \
		float si12 = *sigma12; \
		float ga23 = *gamma23; \
		float si23 = *sigma23; \
		float temp1; \
		float temp2; \
		float temp3; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
\
		alpha1 += 1; \
		alpha2 += 1; \
		alpha3 += 1; \
	} \
}

#define MAC_Apply_G_mx3_asd( m_A, \
                             gamma12, \
                             sigma12, \
                             gamma23, \
                             sigma23, \
                             a1, inc_a1, \
                             a2, inc_a2, \
                             a3, inc_a3 ) \
{\
	int              n_iter16 = m_A / ( 2 * 8 ); \
	int              n_left16 = m_A % ( 2 * 8 ); \
	int              n_iter2  = n_left16 / ( 2 * 1 ); \
	int              n_left   = n_left16 % ( 2 * 1 ); \
	int              i; \
\
	const int        step_a1 = inc_a1 * 2; \
	const int        step_a2 = inc_a1 * 2; \
	const int        step_a3 = inc_a1 * 2; \
\
	double* restrict alpha1 = a1; \
	double* restrict alpha2 = a2; \
	double* restrict alpha3 = a3; \
\
	v2df_t           a1v, a2v, a3v; \
	v2df_t           g12v, s12v; \
	v2df_t           g23v, s23v; \
	v2df_t           t1v, t2v; \
\
	g12v.v = _mm_loaddup_pd( gamma12 ); \
	s12v.v = _mm_loaddup_pd( sigma12 ); \
	g23v.v = _mm_loaddup_pd( gamma23 ); \
	s23v.v = _mm_loaddup_pd( sigma23 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		double ga12 = *gamma12; \
		double si12 = *sigma12; \
		double ga23 = *gamma23; \
		double si23 = *sigma23; \
		double temp1; \
		double temp2; \
		double temp3; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
	} \
}

#define MAC_Apply_G_mx3_asc( m_A, \
                             gamma12, \
                             sigma12, \
                             gamma23, \
                             sigma23, \
                             a1, inc_a1, \
                             a2, inc_a2, \
                             a3, inc_a3 ) \
{ \
	int                n_iter16 = m_A / ( 2 * 8 ); \
	int                n_left16 = m_A % ( 2 * 8 ); \
	int                n_iter2  = n_left16 / ( 2 * 1 ); \
	int                n_left   = n_left16 % ( 2 * 1 ); \
	int                i; \
\
	const int          step_a1 = inc_a1 * 2; \
	const int          step_a2 = inc_a1 * 2; \
	const int          step_a3 = inc_a1 * 2; \
\
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex* restrict alpha3 = a3; \
\
	v4sf_t             a1v, a2v, a3v; \
	v4sf_t             g12v, s12v; \
	v4sf_t             g23v, s23v; \
	v4sf_t             t1v, t2v; \
\
	g12v.v = _mm_load1_ps( gamma12 ); \
	s12v.v = _mm_load1_ps( sigma12 ); \
	g23v.v = _mm_load1_ps( gamma23 ); \
	s23v.v = _mm_load1_ps( sigma23 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		float ga12 = *gamma12; \
		float si12 = *sigma12; \
		float ga23 = *gamma23; \
		float si23 = *sigma23; \
		scomplex temp1; \
		scomplex temp2; \
		scomplex temp3; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real = temp1.real * ga12 + temp2.real * si12; \
		alpha2->real = temp2.real * ga12 - temp1.real * si12; \
\
		alpha1->imag = temp1.imag * ga12 + temp2.imag * si12; \
		alpha2->imag = temp2.imag * ga12 - temp1.imag * si12; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23 + temp3.real * si23; \
		alpha3->real = temp3.real * ga23 - temp2.real * si23; \
\
		alpha2->imag = temp2.imag * ga23 + temp3.imag * si23; \
		alpha3->imag = temp3.imag * ga23 - temp2.imag * si23; \
	} \
}

#define MAC_Apply_G_mx3_asz( m_A, \
                             gamma12, \
                             sigma12, \
                             gamma23, \
                             sigma23, \
                             a1, inc_a1, \
                             a2, inc_a2, \
                             a3, inc_a3 ) \
{\
	int                n_iter = m_A / 8; \
	int                n_left = m_A % 8; \
	int                i; \
\
	const int          step_a1 = inc_a1 * 1; \
	const int          step_a2 = inc_a1 * 1; \
	const int          step_a3 = inc_a1 * 1; \
\
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex* restrict alpha3 = a3; \
\
	v2df_t             a1v, a2v, a3v; \
	v2df_t             g12v, s12v; \
	v2df_t             g23v, s23v; \
	v2df_t             t1v, t2v; \
\
	g12v.v = _mm_loaddup_pd( gamma12 ); \
	s12v.v = _mm_loaddup_pd( sigma12 ); \
	g23v.v = _mm_loaddup_pd( gamma23 ); \
	s23v.v = _mm_loaddup_pd( sigma23 ); \
\
	for ( i = 0; i < n_iter; ++i ) \
	{ \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha2 += step_a2; \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12v.v + a2v.v * s12v.v; \
		a2v.v = a2v.v * g12v.v - t1v.v * s12v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23v.v + a3v.v * s23v.v; \
		a3v.v = a3v.v * g23v.v - t2v.v * s23v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
}

#endif
