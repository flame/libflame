/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Apply_GTG( FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon1, FLA_Obj delta2 );
FLA_Error FLA_Apply_GTG_ops( float*  gamma,
                             float*  sigma,
                             float*  delta1,
                             float*  epsilon1,
                             float*  delta2 );
FLA_Error FLA_Apply_GTG_opd( double* gamma,
                             double* sigma,
                             double* delta1,
                             double* epsilon1,
                             double* delta2 );

#define MAC_Apply_GTG_ops( gamma, sigma, delta1, epsilon, delta2 ) \
{ \
	float  g, s; \
	float  d1, e, d2; \
	float  g2, s2, tgse; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e  = *(epsilon); \
	d2 = *(delta2); \
\
	g2 = g * g; \
	s2 = s * s; \
	tgse = 2.0 * g * s * e; \
\
	*(delta1)  = g2 * d1 + tgse + s2 * d2; \
	*(delta2)  = s2 * d1 - tgse + g2 * d2; \
	*(epsilon) = g * s * (d2 - d1) + e * (g2 - s2); \
}

#define MAC_Apply_GTG_opd( gamma, sigma, delta1, epsilon, delta2 ) \
{ \
/*
	double g, s; \
	double d1, e, d2; \
	double t, st; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e  = *(epsilon); \
	d2 = *(delta2); \
\
	t   = s * ( d2 - d1 ) + 2.0 * g * e; \
	st  = s * t; \
	e   = g * t - e; \
	d1  = st + d1; \
	d2  = d2 - st; \
\
	*(delta1)  = d1; \
	*(epsilon) = e; \
	*(delta2)  = d2; \
*/ \
	double g, s; \
	double d1, e, d2; \
	double g2, s2, tgse; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e  = *(epsilon); \
	d2 = *(delta2); \
\
	g2 = g * g; \
	s2 = s * s; \
	tgse = 2.0 * g * s * e; \
\
	*(delta1)  = g2 * d1 + tgse + s2 * d2; \
	*(delta2)  = s2 * d1 - tgse + g2 * d2; \
	*(epsilon) = g * s * (d2 - d1) + e * (g2 - s2); \
\
/*
	double g, s; \
	double d1, e, d2; \
	double g2, s2; \
	double st; \
\
	g = *(gamma); \
	s = *(sigma); \
\
	d1 = *(delta1); \
	e  = *(epsilon); \
	d2 = *(delta2); \
\
	g2 = g * g; \
	s2 = s * s; \
	st = s2 * (d2 - d1) + 2.0 * g * s * e; \
\
	*(delta1)  = st + d1; \
	*(delta2)  = d2 - st; \
	*(epsilon) = g * s * (d2 - d1) + e * (g2 - s2); \
*/ \
}

