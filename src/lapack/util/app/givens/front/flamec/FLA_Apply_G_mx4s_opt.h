/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_G_mx4s_ops( m_A, \
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
{ \
	float              ga23_k1 = *gamma23_k1; \
	float              si23_k1 = *sigma23_k1; \
	float              ga34_k1 = *gamma34_k1; \
	float              si34_k1 = *sigma34_k1; \
	float              ga12_k2 = *gamma12_k2; \
	float              si12_k2 = *sigma12_k2; \
	float              ga23_k2 = *gamma23_k2; \
	float              si23_k2 = *sigma23_k2; \
	float*    restrict alpha1 = a1; \
	float*    restrict alpha2 = a2; \
	float*    restrict alpha3 = a3; \
	float*    restrict alpha4 = a4; \
	float              temp1; \
	float              temp2; \
	float              temp3; \
	float              temp4; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
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
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
		alpha4 += inc_a4; \
	} \
}

#define MAC_Apply_G_mx4s_opc( m_A, \
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
{ \
	float              ga23_k1 = *gamma23_k1; \
	float              si23_k1 = *sigma23_k1; \
	float              ga34_k1 = *gamma34_k1; \
	float              si34_k1 = *sigma34_k1; \
	float              ga12_k2 = *gamma12_k2; \
	float              si12_k2 = *sigma12_k2; \
	float              ga23_k2 = *gamma23_k2; \
	float              si23_k2 = *sigma23_k2; \
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex* restrict alpha3 = a3; \
	scomplex* restrict alpha4 = a4; \
	scomplex           temp1; \
	scomplex           temp2; \
	scomplex           temp3; \
	scomplex           temp4; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k1 + temp3.real * si23_k1; \
		alpha2->imag = temp2.imag * ga23_k1 + temp3.imag * si23_k1; \
\
		alpha3->real = temp3.real * ga23_k1 - temp2.real * si23_k1; \
		alpha3->imag = temp3.imag * ga23_k1 - temp2.imag * si23_k1; \
\
		temp3 = *alpha3; \
		temp4 = *alpha4; \
\
		alpha3->real = temp3.real * ga34_k1 + temp4.real * si34_k1; \
		alpha3->imag = temp3.imag * ga34_k1 + temp4.imag * si34_k1; \
\
		alpha4->real = temp4.real * ga34_k1 - temp3.real * si34_k1; \
		alpha4->imag = temp4.imag * ga34_k1 - temp3.imag * si34_k1; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real = temp1.real * ga12_k2 + temp2.real * si12_k2; \
		alpha1->imag = temp1.imag * ga12_k2 + temp2.imag * si12_k2; \
\
		alpha2->real = temp2.real * ga12_k2 - temp1.real * si12_k2; \
		alpha2->imag = temp2.imag * ga12_k2 - temp1.imag * si12_k2; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k2 + temp3.real * si23_k2; \
		alpha2->imag = temp2.imag * ga23_k2 + temp3.imag * si23_k2; \
\
		alpha3->real = temp3.real * ga23_k2 - temp2.real * si23_k2; \
		alpha3->imag = temp3.imag * ga23_k2 - temp2.imag * si23_k2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
		alpha4 += inc_a4; \
	} \
}

#define MAC_Apply_G_mx4s_opd( m_A, \
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
{ \
	double             ga23_k1 = *gamma23_k1; \
	double             si23_k1 = *sigma23_k1; \
	double             ga34_k1 = *gamma34_k1; \
	double             si34_k1 = *sigma34_k1; \
	double             ga12_k2 = *gamma12_k2; \
	double             si12_k2 = *sigma12_k2; \
	double             ga23_k2 = *gamma23_k2; \
	double             si23_k2 = *sigma23_k2; \
	double*   restrict alpha1 = a1; \
	double*   restrict alpha2 = a2; \
	double*   restrict alpha3 = a3; \
	double*   restrict alpha4 = a4; \
	double             temp1; \
	double             temp2; \
	double             temp3; \
	double             temp4; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
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
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
		alpha4 += inc_a4; \
	} \
}

#define MAC_Apply_G_mx4s_opz( m_A, \
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
{ \
	double             ga23_k1 = *gamma23_k1; \
	double             si23_k1 = *sigma23_k1; \
	double             ga34_k1 = *gamma34_k1; \
	double             si34_k1 = *sigma34_k1; \
	double             ga12_k2 = *gamma12_k2; \
	double             si12_k2 = *sigma12_k2; \
	double             ga23_k2 = *gamma23_k2; \
	double             si23_k2 = *sigma23_k2; \
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex* restrict alpha3 = a3; \
	dcomplex* restrict alpha4 = a4; \
	dcomplex           temp1; \
	dcomplex           temp2; \
	dcomplex           temp3; \
	dcomplex           temp4; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k1 + temp3.real * si23_k1; \
		alpha2->imag = temp2.imag * ga23_k1 + temp3.imag * si23_k1; \
\
		alpha3->real = temp3.real * ga23_k1 - temp2.real * si23_k1; \
		alpha3->imag = temp3.imag * ga23_k1 - temp2.imag * si23_k1; \
\
		temp3 = *alpha3; \
		temp4 = *alpha4; \
\
		alpha3->real = temp3.real * ga34_k1 + temp4.real * si34_k1; \
		alpha3->imag = temp3.imag * ga34_k1 + temp4.imag * si34_k1; \
\
		alpha4->real = temp4.real * ga34_k1 - temp3.real * si34_k1; \
		alpha4->imag = temp4.imag * ga34_k1 - temp3.imag * si34_k1; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real = temp1.real * ga12_k2 + temp2.real * si12_k2; \
		alpha1->imag = temp1.imag * ga12_k2 + temp2.imag * si12_k2; \
\
		alpha2->real = temp2.real * ga12_k2 - temp1.real * si12_k2; \
		alpha2->imag = temp2.imag * ga12_k2 - temp1.imag * si12_k2; \
\
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real = temp2.real * ga23_k2 + temp3.real * si23_k2; \
		alpha2->imag = temp2.imag * ga23_k2 + temp3.imag * si23_k2; \
\
		alpha3->real = temp3.real * ga23_k2 - temp2.real * si23_k2; \
		alpha3->imag = temp3.imag * ga23_k2 - temp2.imag * si23_k2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
		alpha4 += inc_a4; \
	} \
}

