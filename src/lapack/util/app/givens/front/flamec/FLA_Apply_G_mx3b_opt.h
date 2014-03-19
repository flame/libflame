/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_G_mx3b_ops( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{ \
	float              ga12   = *gamma12; \
	float              si12   = *sigma12; \
	float              ga23   = *gamma23; \
	float              si23   = *sigma23; \
	float*    restrict alpha1 = a1; \
	float*    restrict alpha2 = a2; \
	float*    restrict alpha3 = a3; \
	float              temp1; \
	float              temp2; \
	float              temp3; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
	} \
}

#define MAC_Apply_G_mx3b_opc( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{ \
	float              ga12   = *gamma12; \
	float              si12   = *sigma12; \
	float              ga23   = *gamma23; \
	float              si23   = *sigma23; \
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex* restrict alpha3 = a3; \
	scomplex           temp1; \
	scomplex           temp2; \
	scomplex           temp3; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real =  ga23 * temp2.real + si23 * temp3.real; \
		alpha2->imag =  ga23 * temp2.imag + si23 * temp3.imag; \
\
		alpha3->real = -si23 * temp2.real + ga23 * temp3.real; \
		alpha3->imag = -si23 * temp2.imag + ga23 * temp3.imag; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real =  ga12 * temp1.real + si12 * temp2.real; \
		alpha1->imag =  ga12 * temp1.imag + si12 * temp2.imag; \
\
		alpha2->real = -si12 * temp1.real + ga12 * temp2.real; \
		alpha2->imag = -si12 * temp1.imag + ga12 * temp2.imag; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
	} \
}

#define MAC_Apply_G_mx3b_opd( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{ \
	double             ga12   = *gamma12; \
	double             si12   = *sigma12; \
	double             ga23   = *gamma23; \
	double             si23   = *sigma23; \
	double*   restrict alpha1 = a1; \
	double*   restrict alpha2 = a2; \
	double*   restrict alpha3 = a3; \
	double             temp1; \
	double             temp2; \
	double             temp3; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		*alpha2 = temp2 * ga23 + temp3 * si23; \
		*alpha3 = temp3 * ga23 - temp2 * si23; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 = temp1 * ga12 + temp2 * si12; \
		*alpha2 = temp2 * ga12 - temp1 * si12; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
	} \
}

#define MAC_Apply_G_mx3b_opz( m_A, \
                              gamma12, \
                              sigma12, \
                              gamma23, \
                              sigma23, \
                              a1, inc_a1, \
                              a2, inc_a2, \
                              a3, inc_a3 ) \
{ \
	double             ga12   = *gamma12; \
	double             si12   = *sigma12; \
	double             ga23   = *gamma23; \
	double             si23   = *sigma23; \
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex* restrict alpha3 = a3; \
	dcomplex           temp1; \
	dcomplex           temp2; \
	dcomplex           temp3; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp2 = *alpha2; \
		temp3 = *alpha3; \
\
		alpha2->real =  ga23 * temp2.real + si23 * temp3.real; \
		alpha2->imag =  ga23 * temp2.imag + si23 * temp3.imag; \
\
		alpha3->real = -si23 * temp2.real + ga23 * temp3.real; \
		alpha3->imag = -si23 * temp2.imag + ga23 * temp3.imag; \
\
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		alpha1->real =  ga12 * temp1.real + si12 * temp2.real; \
		alpha1->imag =  ga12 * temp1.imag + si12 * temp2.imag; \
\
		alpha2->real = -si12 * temp1.real + ga12 * temp2.real; \
		alpha2->imag = -si12 * temp1.imag + ga12 * temp2.imag; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
		alpha3 += inc_a3; \
	} \
}

