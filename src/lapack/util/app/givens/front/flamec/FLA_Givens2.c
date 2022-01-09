/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Givens2( FLA_Obj chi_1, FLA_Obj chi_2, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj chi_1_new )
/*
  Compute a Givens rotation

    G  =  / gamma  -sigma \
          \ sigma   gamma /

  by computing gamma and sigma such that the following is satisfied:
  

    G' / chi_1 \ = / alpha \
       \ chi_2 /   \   0   /

    where

      alpha  = +/- || x ||_2

      x      =  / chi_1 \
                \ chi_2 /

  Note that this is equivalent to computing gamma and sigma such that

    ( chi_1  chi_2 ) G = ( alpha  0 )

    where

      alpha  = +/- || x ||_2

      x      =  ( chi_1  chi_2 )

  Upon completion, chi_1_new is overwritten with alpha.

  -FGVZ
*/
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( chi_1 );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Givens2_check( chi_1, chi_2, gamma, sigma, chi_1_new );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float* chi_1_p     = ( float* ) FLA_FLOAT_PTR( chi_1 );
			float* chi_2_p     = ( float* ) FLA_FLOAT_PTR( chi_2 );
			float* gamma_p     = ( float* ) FLA_FLOAT_PTR( gamma );
			float* sigma_p     = ( float* ) FLA_FLOAT_PTR( sigma );
			float* chi_1_new_p = ( float* ) FLA_FLOAT_PTR( chi_1_new );

			FLA_Givens2_ops( chi_1_p,
			                 chi_2_p,
			                 gamma_p,
			                 sigma_p,
			                 chi_1_new_p );

			break;
		}

		case FLA_DOUBLE:
		{
			double* chi_1_p     = ( double* ) FLA_DOUBLE_PTR( chi_1 );
			double* chi_2_p     = ( double* ) FLA_DOUBLE_PTR( chi_2 );
			double* gamma_p     = ( double* ) FLA_DOUBLE_PTR( gamma );
			double* sigma_p     = ( double* ) FLA_DOUBLE_PTR( sigma );
			double* chi_1_new_p = ( double* ) FLA_DOUBLE_PTR( chi_1_new );

			FLA_Givens2_opd( chi_1_p,
			                 chi_2_p,
			                 gamma_p,
			                 sigma_p,
			                 chi_1_new_p );

			break;
		}

	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Givens2_ops( float*  chi_1,
                           float*  chi_2,
                           float*  gamma,
                           float*  sigma,
                           float*  chi_1_new )
{
	return FLA_SUCCESS;
}

FLA_Error FLA_Givens2_opd( double* chi_1,
                           double* chi_2,
                           double* gamma,
                           double* sigma,
                           double* chi_1_new )
{
	double zero = 0.0;
	double one  = 1.0;
	double two  = 2.0;

	double f = *chi_1;
	double g = *chi_2;
	double cs;
	double sn;
	double r;

	integer    count, i;
	double eps, f1, g1, safmin, safmin2, safmax2, scale;
	double base;

	safmin  = FLA_Mach_params_opd( FLA_MACH_SFMIN );
	eps     = FLA_Mach_params_opd( FLA_MACH_EPS );
	base    = FLA_Mach_params_opd( FLA_MACH_BASE );

	safmin2 = pow( base, ( double )(( integer )( log( safmin / eps ) /
	                                         log( base ) /
	                                         two ) ) );
	safmax2 = one / safmin2;

	if ( g == zero )
	{
		cs = one;
		sn = zero;
		r  = f;
	}
	else if ( f == zero )
	{
		cs = zero;
		sn = one;
		r  = g;
	}
	else
	{
		f1 = f;
		g1 = g;
		scale = max( fabs( f1 ), fabs( g1 ) );

		if ( scale >= safmax2 )
		{
			count = 0;
			do
			{
				++count;
				f1 = f1 * safmin2;
				g1 = g1 * safmin2;
				scale = max( fabs( f1 ), fabs( g1 ) );
			}
			while ( scale >= safmax2 );

			r  = sqrt( f1 * f1 + g1 * g1 );
			cs = f1 / r;
			sn = g1 / r;

			for ( i = 0; i < count; ++i )
				r = r * safmax2;
		}
		else if ( scale <= safmin2 )
		{
			count = 0;
			do
			{
				++count;
				f1 = f1 * safmax2;
				g1 = g1 * safmax2;
				scale = max( fabs( f1 ), fabs( g1 ) );
			}
			while ( scale <= safmin2 );

			r  = sqrt( f1 * f1 + g1 * g1 );
			cs = f1 / r;
			sn = g1 / r;

			for ( i = 0; i < count; ++i )
				r = r * safmin2;
		}
		else
		{
			r  = sqrt( f1 * f1 + g1 * g1 );
			cs = f1 / r;
			sn = g1 / r;
		}

		if ( fabs( f ) > fabs ( g ) && cs < zero )
		{
			cs = -cs;
			sn = -sn;
			r  = -r;
		}
	}

	// Save the output values.
	*gamma     = cs;
	*sigma     = sn;
	*chi_1_new = r;

	return FLA_SUCCESS;
}

