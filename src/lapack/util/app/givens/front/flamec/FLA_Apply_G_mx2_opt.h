/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_G_mx2_ops( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	float             ga     = *gamma12; \
	float             si     = *sigma12; \
	float*  restrict  alpha1 = a1; \
	float*  restrict  alpha2 = a2; \
	float             temp1; \
	float             temp2; \
	int               i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 =  ga * temp1 + si * temp2; \
		*alpha2 = -si * temp1 + ga * temp2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}

#define MAC_Apply_G_mx2_opc( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	float              ga12   = *gamma12; \
	float              si12   = *sigma12; \
	scomplex* restrict alpha1 = a1; \
	scomplex* restrict alpha2 = a2; \
	scomplex           temp1; \
	scomplex           temp2; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
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
	} \
}

#define MAC_Apply_G_mx2_opd( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{ \
	double            ga     = *gamma12; \
	double            si     = *sigma12; \
	double* restrict  alpha1 = a1; \
	double* restrict  alpha2 = a2; \
	double            temp1; \
	double            temp2; \
	int               i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
		temp1 = *alpha1; \
		temp2 = *alpha2; \
\
		*alpha1 =  ga * temp1 + si * temp2; \
		*alpha2 = -si * temp1 + ga * temp2; \
\
		alpha1 += inc_a1; \
		alpha2 += inc_a2; \
	} \
}

#define MAC_Apply_G_mx2_opz( m_A, \
                             gamma12, \
                             sigma12, \
                             a1, inc_a1, \
                             a2, inc_a2 ) \
{\
	double             ga12   = *gamma12; \
	double             si12   = *sigma12; \
	dcomplex* restrict alpha1 = a1; \
	dcomplex* restrict alpha2 = a2; \
	dcomplex           temp1; \
	dcomplex           temp2; \
	int                i; \
\
	for ( i = 0; i < m_A; ++i ) \
	{ \
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
	} \
}

