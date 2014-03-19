/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_GT_2x2_ops( gamma, sigma, epsilon1, delta2, beta, epsilon2 ) \
{ \
	float g, s; \
	float e1, d2, e2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	e1 = *(epsilon1); \
	d2 = *(delta2); \
	e2 = *(epsilon2); \
\
	*(epsilon1)  =  g * e1 + s * d2; \
	*(delta2)    = -s * e1 + g * d2; \
\
	*(beta)      = s * e2; \
	*(epsilon2)  = g * e2; \
}

#define MAC_Apply_GT_2x2_opd( gamma, sigma, epsilon1, delta2, beta, epsilon2 ) \
{ \
	double g, s; \
	double e1, d2, e2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	e1 = *(epsilon1); \
	d2 = *(delta2); \
	e2 = *(epsilon2); \
\
	*(epsilon1)  =  g * e1 + s * d2; \
	*(delta2)    = -s * e1 + g * d2; \
\
	*(beta)      = s * e2; \
	*(epsilon2)  = g * e2; \
}

#define MAC_Apply_GT_2x1_ops( gamma, sigma, epsilon1, delta2 ) \
{ \
	float g, s; \
	float e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(epsilon1)  =  g * e1 + s * d2; \
	*(delta2)    = -s * e1 + g * d2; \
}

#define MAC_Apply_GT_2x1_opd( gamma, sigma, epsilon1, delta2 ) \
{ \
	double g, s; \
	double e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(epsilon1)  =  g * e1 + s * d2; \
	*(delta2)    = -s * e1 + g * d2; \
}

