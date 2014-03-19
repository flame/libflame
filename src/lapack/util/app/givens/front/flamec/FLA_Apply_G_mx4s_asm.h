/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


#if FLA_VECTOR_INTRINSIC_TYPE == FLA_NO_INTRINSICS

#define MAC_Apply_G_mx4s_ass MAC_Apply_G_mx4s_ops
#define MAC_Apply_G_mx4s_asd MAC_Apply_G_mx4s_opd
#define MAC_Apply_G_mx4s_asc MAC_Apply_G_mx4s_opc
#define MAC_Apply_G_mx4s_asz MAC_Apply_G_mx4s_opz

#elif FLA_VECTOR_INTRINSIC_TYPE == FLA_SSE_INTRINSICS

#define MAC_Apply_G_mx4s_ass( m_A, \
                              gamma23_k1, \
                              sigma23_k1, \
                              gamma34_k1, \
                              sigma34_k1, \
                              gamma12_k2, \
                              sigma12_k2, \
                              gamma23_k2, \
                              sigma23_k2, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3, \
                              a4, inc_a4 ) \
{\
	int                n_iter32 = m_A / ( 4 * 8 ); \
	int                n_left32 = m_A % ( 4 * 8 ); \
	int                n_iter4  = n_left32 / ( 4 * 1 ); \
	int                n_left   = n_left32 % ( 4 * 1 ); \
	int                i; \
\
	const int          step_a1 = inc_a1 * 4; \
	const int          step_a2 = inc_a2 * 4; \
	const int          step_a3 = inc_a3 * 4; \
	const int          step_a4 = inc_a4 * 4; \
\
	float*    restrict alpha1 = a1; \
	float*    restrict alpha2 = a2; \
	float*    restrict alpha3 = a3; \
	float*    restrict alpha4 = a4; \
\
	v4sf_t             a1v, a2v, a3v, a4v; \
	v4sf_t             b1v, b2v, b3v, b4v; \
	v4sf_t             g23_k1v, s23_k1v; \
	v4sf_t             g34_k1v, s34_k1v; \
	v4sf_t             g12_k2v, s12_k2v; \
	v4sf_t             g23_k2v, s23_k2v; \
	v4sf_t             t1v, t2v, t3v; \
\
	g23_k1v.v = _mm_load1_ps( gamma23_k1 ); \
	s23_k1v.v = _mm_load1_ps( sigma23_k1 ); \
	g34_k1v.v = _mm_load1_ps( gamma34_k1 ); \
	s34_k1v.v = _mm_load1_ps( sigma34_k1 ); \
	g12_k2v.v = _mm_load1_ps( gamma12_k2 ); \
	s12_k2v.v = _mm_load1_ps( sigma12_k2 ); \
	g23_k2v.v = _mm_load1_ps( gamma23_k2 ); \
	s23_k2v.v = _mm_load1_ps( sigma23_k2 ); \
\
	for ( i = 0; i < n_iter32; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a3) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
\
/* ----------------------------------------------------------- */ \
	} \
\
	for ( i = 0; i < n_iter4; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
		float              ga23_k1 = *gamma23_k1; \
		float              si23_k1 = *sigma23_k1; \
		float              ga34_k1 = *gamma34_k1; \
		float              si34_k1 = *sigma34_k1; \
		float              ga12_k2 = *gamma12_k2; \
		float              si12_k2 = *sigma12_k2; \
		float              ga23_k2 = *gamma23_k2; \
		float              si23_k2 = *sigma23_k2; \
		float              temp1; \
		float              temp2; \
		float              temp3; \
		float              temp4; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23_k1 + temp3 * si23_k1; \
		*alpha3 = temp3 * ga23_k1 - temp2 * si23_k1; \
\
		temp3 = *alpha3; \
		temp4 = *alpha4; \
\
		*alpha3 = temp3 * ga34_k1 + temp4 * si34_k1; \
		*alpha4 = temp4 * ga34_k1 - temp3 * si34_k1; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12_k2 + temp2 * si12_k2; \
		*alpha2 = temp2 * ga12_k2 - temp1 * si12_k2; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23_k2 + temp3 * si23_k2; \
		*alpha3 = temp3 * ga23_k2 - temp2 * si23_k2; \
\
		alpha1 += 1; \
		alpha2 += 1; \
		alpha3 += 1; \
		alpha4 += 1; \
	} \
}

#define MAC_Apply_G_mx4s_asd( m_A, \
                              gamma23_k1, \
                              sigma23_k1, \
                              gamma34_k1, \
                              sigma34_k1, \
                              gamma12_k2, \
                              sigma12_k2, \
                              gamma23_k2, \
                              sigma23_k2, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3, \
                              a4, inc_a4 ) \
{\
	int                n_iter16 = m_A / ( 2 * 8 ); \
	int                n_left16 = m_A % ( 2 * 8 ); \
	int                n_iter2  = n_left16 / ( 2 * 1 ); \
	int                n_left   = n_left16 % ( 2 * 1 ); \
	int                i; \
\
	const int          step_a1 = inc_a1 * 2; \
	const int          step_a2 = inc_a2 * 2; \
	const int          step_a3 = inc_a3 * 2; \
	const int          step_a4 = inc_a4 * 2; \
\
	double*   restrict alpha1 = a1; \
	double*   restrict alpha2 = a2; \
	double*   restrict alpha3 = a3; \
	double*   restrict alpha4 = a4; \
\
	v2df_t             a1v, a2v, a3v, a4v; \
	v2df_t             b1v, b2v, b3v, b4v; \
	v2df_t             g23_k1v, s23_k1v; \
	v2df_t             g34_k1v, s34_k1v; \
	v2df_t             g12_k2v, s12_k2v; \
	v2df_t             g23_k2v, s23_k2v; \
	v2df_t             t1v, t2v, t3v; \
\
	g23_k1v.v = _mm_loaddup_pd( gamma23_k1 ); \
	s23_k1v.v = _mm_loaddup_pd( sigma23_k1 ); \
	g34_k1v.v = _mm_loaddup_pd( gamma34_k1 ); \
	s34_k1v.v = _mm_loaddup_pd( sigma34_k1 ); \
	g12_k2v.v = _mm_loaddup_pd( gamma12_k2 ); \
	s12_k2v.v = _mm_loaddup_pd( sigma12_k2 ); \
	g23_k2v.v = _mm_loaddup_pd( gamma23_k2 ); \
	s23_k2v.v = _mm_loaddup_pd( sigma23_k2 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a3) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
\
/* ----------------------------------------------------------- */ \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		double             ga23_k1 = *gamma23_k1; \
		double             si23_k1 = *sigma23_k1; \
		double             ga34_k1 = *gamma34_k1; \
		double             si34_k1 = *sigma34_k1; \
		double             ga12_k2 = *gamma12_k2; \
		double             si12_k2 = *sigma12_k2; \
		double             ga23_k2 = *gamma23_k2; \
		double             si23_k2 = *sigma23_k2; \
		double             temp1; \
		double             temp2; \
		double             temp3; \
		double             temp4; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23_k1 + temp3 * si23_k1; \
		*alpha3 = temp3 * ga23_k1 - temp2 * si23_k1; \
\
		temp3 = *alpha3; \
		temp4 = *alpha4; \
\
		*alpha3 = temp3 * ga34_k1 + temp4 * si34_k1; \
		*alpha4 = temp4 * ga34_k1 - temp3 * si34_k1; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12_k2 + temp2 * si12_k2; \
		*alpha2 = temp2 * ga12_k2 - temp1 * si12_k2; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23_k2 + temp3 * si23_k2; \
		*alpha3 = temp3 * ga23_k2 - temp2 * si23_k2; \
\
	} \
}

#define MAC_Apply_G_mx4s_asc( m_A, \
                              gamma23_k1, \
                              sigma23_k1, \
                              gamma34_k1, \
                              sigma34_k1, \
                              gamma12_k2, \
                              sigma12_k2, \
                              gamma23_k2, \
                              sigma23_k2, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3, \
                              a4, inc_a4 ) \
{\
	int                n_iter16 = m_A / ( 2 * 8 ); \
	int                n_left16 = m_A % ( 2 * 8 ); \
	int                n_iter2  = n_left16 / ( 2 * 1 ); \
	int                n_left   = n_left16 % ( 2 * 1 ); \
	int                i; \
\
	const int          step_a1 = inc_a1 * 2; \
	const int          step_a2 = inc_a2 * 2; \
	const int          step_a3 = inc_a3 * 2; \
	const int          step_a4 = inc_a4 * 2; \
\
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex* restrict alpha3 = a3; \
	scomplex* restrict alpha4 = a4; \
\
	v4sf_t             a1v, a2v, a3v, a4v; \
	v4sf_t             b1v, b2v, b3v, b4v; \
	v4sf_t             g23_k1v, s23_k1v; \
	v4sf_t             g34_k1v, s34_k1v; \
	v4sf_t             g12_k2v, s12_k2v; \
	v4sf_t             g23_k2v, s23_k2v; \
	v4sf_t             t1v, t2v, t3v; \
\
	g23_k1v.v = _mm_load1_ps( gamma23_k1 ); \
	s23_k1v.v = _mm_load1_ps( sigma23_k1 ); \
	g34_k1v.v = _mm_load1_ps( gamma34_k1 ); \
	s34_k1v.v = _mm_load1_ps( sigma34_k1 ); \
	g12_k2v.v = _mm_load1_ps( gamma12_k2 ); \
	s12_k2v.v = _mm_load1_ps( sigma12_k2 ); \
	g23_k2v.v = _mm_load1_ps( gamma23_k2 ); \
	s23_k2v.v = _mm_load1_ps( sigma23_k2 ); \
\
	for ( i = 0; i < n_iter16; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a3) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_ps( ( float* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_ps( ( float* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
\
		_mm_store_ps( ( float* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
\
/* ----------------------------------------------------------- */ \
	} \
\
	for ( i = 0; i < n_iter2; ++i ) \
	{ \
\
		a2v.v = _mm_load_ps( ( float* )alpha2 ); \
		a3v.v = _mm_load_ps( ( float* )alpha3 ); \
		a4v.v = _mm_load_ps( ( float* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_ps( ( float* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_ps( ( float* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_ps( ( float* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_ps( ( float* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_ps( ( float* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
\
	if ( n_left == 1 ) \
	{ \
		float             ga23_k1 = *gamma23_k1; \
		float             si23_k1 = *sigma23_k1; \
		float             ga34_k1 = *gamma34_k1; \
		float             si34_k1 = *sigma34_k1; \
		float             ga12_k2 = *gamma12_k2; \
		float             si12_k2 = *sigma12_k2; \
		float             ga23_k2 = *gamma23_k2; \
		float             si23_k2 = *sigma23_k2; \
		scomplex          temp1; \
		scomplex          temp2; \
		scomplex          temp3; \
		scomplex          temp4; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k1 + temp3.real * si23_k1; \
		alpha3->real = temp3.real * ga23_k1 - temp2.real * si23_k1; \
\
		alpha2->imag = temp2.imag * ga23_k1 + temp3.imag * si23_k1; \
		alpha3->imag = temp3.imag * ga23_k1 - temp2.imag * si23_k1; \
\
		temp3 = *alpha3; \
		temp4 = *alpha4; \
\
		alpha3->real = temp3.real * ga34_k1 + temp4.real * si34_k1; \
		alpha4->real = temp4.real * ga34_k1 - temp3.real * si34_k1; \
\
		alpha3->imag = temp3.imag * ga34_k1 + temp4.imag * si34_k1; \
		alpha4->imag = temp4.imag * ga34_k1 - temp3.imag * si34_k1; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real = temp1.real * ga12_k2 + temp2.real * si12_k2; \
		alpha2->real = temp2.real * ga12_k2 - temp1.real * si12_k2; \
\
		alpha1->imag = temp1.imag * ga12_k2 + temp2.imag * si12_k2; \
		alpha2->imag = temp2.imag * ga12_k2 - temp1.imag * si12_k2; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k2 + temp3.real * si23_k2; \
		alpha3->real = temp3.real * ga23_k2 - temp2.real * si23_k2; \
\
		alpha2->imag = temp2.imag * ga23_k2 + temp3.imag * si23_k2; \
		alpha3->imag = temp3.imag * ga23_k2 - temp2.imag * si23_k2; \
\
	} \
}

#define MAC_Apply_G_mx4s_asz( m_A, \
                              gamma23_k1, \
                              sigma23_k1, \
                              gamma34_k1, \
                              sigma34_k1, \
                              gamma12_k2, \
                              sigma12_k2, \
                              gamma23_k2, \
                              sigma23_k2, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3, \
                              a4, inc_a4 ) \
{\
	int                n_iter = m_A / 8; \
	int                n_left = m_A % 8; \
	int                i; \
\
	const int          step_a1 = inc_a1 * 1; \
	const int          step_a2 = inc_a2 * 1; \
	const int          step_a3 = inc_a3 * 1; \
	const int          step_a4 = inc_a4 * 1; \
\
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex* restrict alpha3 = a3; \
	dcomplex* restrict alpha4 = a4; \
\
	v2df_t             a1v, a2v, a3v, a4v; \
	v2df_t             b1v, b2v, b3v, b4v; \
	v2df_t             g23_k1v, s23_k1v; \
	v2df_t             g34_k1v, s34_k1v; \
	v2df_t             g12_k2v, s12_k2v; \
	v2df_t             g23_k2v, s23_k2v; \
	v2df_t             t1v, t2v, t3v; \
\
	g23_k1v.v = _mm_loaddup_pd( gamma23_k1 ); \
	s23_k1v.v = _mm_loaddup_pd( sigma23_k1 ); \
	g34_k1v.v = _mm_loaddup_pd( gamma34_k1 ); \
	s34_k1v.v = _mm_loaddup_pd( sigma34_k1 ); \
	g12_k2v.v = _mm_loaddup_pd( gamma12_k2 ); \
	s12_k2v.v = _mm_loaddup_pd( sigma12_k2 ); \
	g23_k2v.v = _mm_loaddup_pd( gamma23_k2 ); \
	s23_k2v.v = _mm_loaddup_pd( sigma23_k2 ); \
\
	for ( i = 0; i < n_iter; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a3) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
		a2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
		a3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
		b2v.v = _mm_load_pd( ( double* )(alpha2 + step_a2) ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		b3v.v = _mm_load_pd( ( double* )(alpha3 + step_a3) ); \
\
/* ----------------------------------------------------------- */ \
\
		b4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k1v.v + b3v.v * s23_k1v.v; \
		b3v.v = b3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
		b1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = b3v.v; \
		b3v.v = t3v.v * g34_k1v.v + b4v.v * s34_k1v.v; \
		b4v.v = b4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, b4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = b1v.v; \
		b1v.v = t1v.v * g12_k2v.v + b2v.v * s12_k2v.v; \
		b2v.v = b2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, b1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = b2v.v; \
		b2v.v = t2v.v * g23_k2v.v + b3v.v * s23_k2v.v; \
		b3v.v = b3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, b2v.v ); \
		alpha2 += step_a2; \
\
		_mm_store_pd( ( double* )alpha3, b3v.v ); \
		alpha3 += step_a3; \
\
/* ----------------------------------------------------------- */ \
	} \
\
	for ( i = 0; i < n_left; ++i ) \
	{ \
\
		a2v.v = _mm_load_pd( ( double* )alpha2 ); \
		a3v.v = _mm_load_pd( ( double* )alpha3 ); \
		a4v.v = _mm_load_pd( ( double* )alpha4 ); \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k1v.v + a3v.v * s23_k1v.v; \
		a3v.v = a3v.v * g23_k1v.v - t2v.v * s23_k1v.v; \
\
		a1v.v = _mm_load_pd( ( double* )alpha1 ); \
\
		t3v.v = a3v.v; \
		a3v.v = t3v.v * g34_k1v.v + a4v.v * s34_k1v.v; \
		a4v.v = a4v.v * g34_k1v.v - t3v.v * s34_k1v.v; \
\
		_mm_store_pd( ( double* )alpha4, a4v.v ); \
		alpha4 += step_a4; \
\
		t1v.v = a1v.v; \
		a1v.v = t1v.v * g12_k2v.v + a2v.v * s12_k2v.v; \
		a2v.v = a2v.v * g12_k2v.v - t1v.v * s12_k2v.v; \
\
		_mm_store_pd( ( double* )alpha1, a1v.v ); \
		alpha1 += step_a1; \
\
		t2v.v = a2v.v; \
		a2v.v = t2v.v * g23_k2v.v + a3v.v * s23_k2v.v; \
		a3v.v = a3v.v * g23_k2v.v - t2v.v * s23_k2v.v; \
\
		_mm_store_pd( ( double* )alpha2, a2v.v ); \
		alpha2 += step_a2; \
		_mm_store_pd( ( double* )alpha3, a3v.v ); \
		alpha3 += step_a3; \
	} \
}

#endif
