/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hevv_2x2( FLA_Obj alpha11, FLA_Obj alpha21, FLA_Obj alpha22,
                        FLA_Obj lambda1, FLA_Obj lambda2,
                        FLA_Obj gamma1,  FLA_Obj sigma1 )
/*
  Compute the eigenvalue decomposition of a 2x2 symmetric matrix A such
  that

    / alpha11 alpha21 \
    \ alpha21 alpha22 /

  is equal to

    / gamma1 -sigma1 \ / lambda1    0    \ / gamma1 -sigma1 \'
    \ sigma1  gamma1 / \    0    lambda2 / \ sigma1  gamma1 /

  Upon completion, lambda1 and lambda2 are overwritten with the
  eigenvalues of larger and smaller absolute values, respectively,
  while gamma1 and sigma1 determine the corresponding pair of 2x1
  eigenvectors.

  This routine is a nearly-verbatim translation of slaev2() and dlaev2()
  from the netlib distribution of LAPACK.

  -FGVZ
*/
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( alpha11 );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  buff_alpha11 = FLA_FLOAT_PTR( alpha11 );
			float*  buff_alpha21 = FLA_FLOAT_PTR( alpha21 );
			float*  buff_alpha22 = FLA_FLOAT_PTR( alpha22 );
			float*  buff_lambda1 = FLA_FLOAT_PTR( lambda1 );
			float*  buff_lambda2 = FLA_FLOAT_PTR( lambda2 );
			float*  buff_gamma1  = FLA_FLOAT_PTR( gamma1 );
			float*  buff_sigma1  = FLA_FLOAT_PTR( sigma1 );

			FLA_Hevv_2x2_ops( buff_alpha11,
			                  buff_alpha21,
			                  buff_alpha22,
			                  buff_lambda1,
			                  buff_lambda2,
			                  buff_gamma1,
			                  buff_sigma1 );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_alpha11 = FLA_DOUBLE_PTR( alpha11 );
			double* buff_alpha21 = FLA_DOUBLE_PTR( alpha21 );
			double* buff_alpha22 = FLA_DOUBLE_PTR( alpha22 );
			double* buff_lambda1 = FLA_DOUBLE_PTR( lambda1 );
			double* buff_lambda2 = FLA_DOUBLE_PTR( lambda2 );
			double* buff_gamma1  = FLA_DOUBLE_PTR( gamma1 );
			double* buff_sigma1  = FLA_DOUBLE_PTR( sigma1 );

			FLA_Hevv_2x2_opd( buff_alpha11,
			                  buff_alpha21,
			                  buff_alpha22,
			                  buff_lambda1,
			                  buff_lambda2,
			                  buff_gamma1,
			                  buff_sigma1 );

			break;
		}

		case FLA_COMPLEX:
		{
			scomplex* buff_alpha11 = FLA_COMPLEX_PTR( alpha11 );
			scomplex* buff_alpha21 = FLA_COMPLEX_PTR( alpha21 );
			scomplex* buff_alpha22 = FLA_COMPLEX_PTR( alpha22 );
			float*    buff_lambda1 = FLA_FLOAT_PTR( lambda1 );
			float*    buff_lambda2 = FLA_FLOAT_PTR( lambda2 );
			float*    buff_gamma1  = FLA_FLOAT_PTR( gamma1 );
			scomplex* buff_sigma1  = FLA_COMPLEX_PTR( sigma1 );

			FLA_Hevv_2x2_opc( buff_alpha11,
			                  buff_alpha21,
			                  buff_alpha22,
			                  buff_lambda1,
			                  buff_lambda2,
			                  buff_gamma1,
			                  buff_sigma1 );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			dcomplex* buff_alpha11 = FLA_DOUBLE_COMPLEX_PTR( alpha11 );
			dcomplex* buff_alpha21 = FLA_DOUBLE_COMPLEX_PTR( alpha21 );
			dcomplex* buff_alpha22 = FLA_DOUBLE_COMPLEX_PTR( alpha22 );
			double*   buff_lambda1 = FLA_DOUBLE_PTR( lambda1 );
			double*   buff_lambda2 = FLA_DOUBLE_PTR( lambda2 );
			double*   buff_gamma1  = FLA_DOUBLE_PTR( gamma1 );
			dcomplex* buff_sigma1  = FLA_DOUBLE_COMPLEX_PTR( sigma1 );

			FLA_Hevv_2x2_opz( buff_alpha11,
			                  buff_alpha21,
			                  buff_alpha22,
			                  buff_lambda1,
			                  buff_lambda2,
			                  buff_gamma1,
			                  buff_sigma1 );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Hevv_2x2_ops( float*    alpha11,
                            float*    alpha21,
                            float*    alpha22,
                            float*    lambda1,
                            float*    lambda2,
                            float*    gamma1,
                            float*    sigma1 )
{
	float  a11, a21, a22;
	float  l1, l2;
	float  g1, s1;
	float  ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn;
	integer    sgn1, sgn2;

	a11 = *alpha11;
	a21 = *alpha21;
	a22 = *alpha22;

	// Compute the eigenvalues.

	sm  = a11 + a22;
	df  = a11 - a22;
	adf = fabs( df );
	tb  = a21 + a21;
	ab  = fabs( tb );

	if ( fabs( a11 ) > fabs( a22 ) )
	{
		acmx = a11;
		acmn = a22;
	}
	else
	{
		acmx = a22;
		acmn = a11;
	}

	if      ( adf > ab ) rt = adf * sqrt( 1.0F + ( ab / adf ) * ( ab / adf ) );
	else if ( adf < ab ) rt = ab  * sqrt( 1.0F + ( adf / ab ) * ( adf / ab ) );
	else                 rt = ab  * sqrt( 2.0F );

	if      ( sm < 0.0F )
	{
		l1 = 0.5F * ( sm - rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
		sgn1 = -1;
	}
	else if ( sm > 0.0F )
	{
		l1 = 0.5F * ( sm + rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
		sgn1 = 1;
	}
	else
	{
		l1 =  0.5F * rt;
		l2 = -0.5F * rt;
		sgn1 = 1;
	}

	*lambda1 = l1;
	*lambda2 = l2;

	// Compute the eigenvector.

	if ( df >= 0.0F )
	{
		cs = df + rt;
		sgn2 = 1;
	}
	else
	{
		cs = df - rt;
		sgn2 = -1;
	}

	acs = fabs( cs );

	if ( acs > ab )
	{
		ct = -tb / cs;
		s1 = 1.0F / sqrt( 1.0F + ct*ct );
		g1 = ct * s1;
	}
	else
	{
		if ( ab == 0.0F )
		{
			g1 = 1.0F;
			s1 = 0.0F;
		}
		else
		{
			tn = -cs / tb;
			g1 = 1.0F / sqrt( 1.0F + tn*tn );
			s1 = tn * g1;
		}
	}
	
	if ( sgn1 == sgn2 )
	{
		tn = g1;
		g1 = -s1;
		s1 = tn;
	}

	*gamma1 = g1;
	*sigma1 = s1;

	return FLA_SUCCESS;
}



FLA_Error FLA_Hevv_2x2_opd( double*   alpha11,
                            double*   alpha21,
                            double*   alpha22,
                            double*   lambda1,
                            double*   lambda2,
                            double*   gamma1,
                            double*   sigma1 )
{
	double a11, a21, a22;
	double l1, l2;
	double g1, s1;
	double ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm, tb, tn;
	integer    sgn1, sgn2;

	a11 = *alpha11;
	a21 = *alpha21;
	a22 = *alpha22;

	// Compute the eigenvalues.

	sm  = a11 + a22;
	df  = a11 - a22;
	adf = fabs( df );
	tb  = a21 + a21;
	ab  = fabs( tb );

	if ( fabs( a11 ) > fabs( a22 ) )
	{
		acmx = a11;
		acmn = a22;
	}
	else
	{
		acmx = a22;
		acmn = a11;
	}

	if      ( adf > ab ) rt = adf * sqrt( 1.0 + pow( ( ab / adf ), 2.0 ) );
	else if ( adf < ab ) rt = ab  * sqrt( 1.0 + pow( ( adf / ab ), 2.0 ) );
	else                 rt = ab  * sqrt( 2.0 );

	if      ( sm < 0.0 )
	{
		l1 = 0.5 * ( sm - rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
		sgn1 = -1;
	}
	else if ( sm > 0.0 )
	{
		l1 = 0.5 * ( sm + rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
		sgn1 = 1;
	}
	else
	{
		l1 =  0.5 * rt;
		l2 = -0.5 * rt;
		sgn1 = 1;
	}

	*lambda1 = l1;
	*lambda2 = l2;

	// Compute the eigenvector.

	if ( df >= 0.0 )
	{
		cs = df + rt;
		sgn2 = 1;
	}
	else
	{
		cs = df - rt;
		sgn2 = -1;
	}

	acs = fabs( cs );

	if ( acs > ab )
	{
		ct = -tb / cs;
		s1 = 1.0 / sqrt( 1.0 + ct*ct );
		g1 = ct * s1;
	}
	else
	{
		if ( ab == 0.0 )
		{
			g1 = 1.0;
			s1 = 0.0;
		}
		else
		{
			tn = -cs / tb;
			g1 = 1.0 / sqrt( 1.0 + tn*tn );
			s1 = tn * g1;
		}
	}
	
	if ( sgn1 == sgn2 )
	{
		tn = g1;
		g1 = -s1;
		s1 = tn;
	}

	*gamma1 = g1;
	*sigma1 = s1;

	return FLA_SUCCESS;
}



FLA_Error FLA_Hevv_2x2_opc( scomplex* alpha11,
                            scomplex* alpha21,
                            scomplex* alpha22,
                            float*    lambda1,
                            float*    lambda2,
                            float*    gamma1,
                            scomplex* sigma1 )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}



FLA_Error FLA_Hevv_2x2_opz( dcomplex* alpha11,
                            dcomplex* alpha21,
                            dcomplex* alpha22,
                            double*   lambda1,
                            double*   lambda2,
                            double*   gamma1,
                            dcomplex* sigma1 )
{
	FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	return FLA_SUCCESS;
}
