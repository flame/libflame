/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#define MAC_Apply_G_2x2_ops( gamma, sigma, delta1, beta, epsilon1, delta2 ) \
{ \
	float g, s; \
	float d1, e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(delta1)    =  g * d1 + s * e1; \
	*(epsilon1)  = -s * d1 + g * e1; \
\
	*(beta)      = s * d2; \
	*(delta2)    = g * d2; \
}

#define MAC_Apply_G_2x2_opd( gamma, sigma, delta1, beta, epsilon1, delta2 ) \
{ \
	double g, s; \
	double d1, e1, d2; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e1 = *(epsilon1); \
	d2 = *(delta2); \
\
	*(delta1)    =  g * d1 + s * e1; \
	*(epsilon1)  = -s * d1 + g * e1; \
\
	*(beta)      = s * d2; \
	*(delta2)    = g * d2; \
}

