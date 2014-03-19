/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Givens2( FLA_Obj chi_1, FLA_Obj chi_2, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj chi_1_new );
FLA_Error FLA_Givens2_ops( float*  chi_1,
                           float*  chi_2,
                           float*  gamma,
                           float*  sigma,
                           float*  chi_1_new );
FLA_Error FLA_Givens2_opd( double* chi_1,
                           double* chi_2,
                           double* gamma,
                           double* sigma,
                           double* chi_1_new );
#define MAC_Givens2_ops( chi_1, chi_2, gamma, sigma, chi_1_new ) \
{ \
	float  chi_1_orig = *(chi_1); \
	float  chi_2_orig = *(chi_2); \
	float  g, s; \
	float  norm_x; \
\
	norm_x = ( float  ) sqrt( ( float ) ( chi_1_orig * chi_1_orig + \
	                                      chi_2_orig * chi_2_orig ) ); \
\
	g = chi_1_orig / norm_x; \
	s = chi_2_orig / norm_x; \
\
	if ( fabs( chi_1_orig ) > fabs( chi_2_orig ) && g < 0.0F ) \
	{ \
		g      = -g; \
		s      = -s; \
		norm_x = -norm_x; \
	} \
\
	*(gamma)     = g; \
	*(sigma)     = s; \
	*(chi_1_new) = norm_x; \
\
}

#define MAC_Givens2_opd( chi_1, chi_2, gamma, sigma, chi_1_new ) \
{ \
	double chi_1_orig = *(chi_1); \
	double chi_2_orig = *(chi_2); \
	double g, s; \
	double norm_x; \
\
	norm_x = ( double ) sqrt( chi_1_orig * chi_1_orig + \
	                          chi_2_orig * chi_2_orig ); \
\
	g = chi_1_orig / norm_x; \
	s = chi_2_orig / norm_x; \
\
	if ( fabs( chi_1_orig ) > fabs( chi_2_orig ) && g < 0.0 ) \
	{ \
		g      = -g; \
		s      = -s; \
		norm_x = -norm_x; \
	} \
\
	*(gamma)     = g; \
	*(sigma)     = s; \
	*(chi_1_new) = norm_x; \
\
}

