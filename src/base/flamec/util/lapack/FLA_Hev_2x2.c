
#include "FLAME.h"

FLA_Error FLA_Hev_2x2( FLA_Obj alpha11, FLA_Obj alpha21, FLA_Obj alpha22,
                       FLA_Obj lambda1, FLA_Obj lambda2 )
/*
  Compute the eigenvalues of a 2x2 symmetric matrix A:

    / alpha11 alpha21 \
    \ alpha21 alpha22 /

  Upon completion, lambda1 and lambda2 are overwritten with the
  eigenvalues of larger and smaller absolute values, respectively.

  This routine is a nearly-verbatim translation of slae2() and dlae2()
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

			FLA_Hev_2x2_ops( buff_alpha11,
			                 buff_alpha21,
			                 buff_alpha22,
			                 buff_lambda1,
			                 buff_lambda2 );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_alpha11 = FLA_DOUBLE_PTR( alpha11 );
			double* buff_alpha21 = FLA_DOUBLE_PTR( alpha21 );
			double* buff_alpha22 = FLA_DOUBLE_PTR( alpha22 );
			double* buff_lambda1 = FLA_DOUBLE_PTR( lambda1 );
			double* buff_lambda2 = FLA_DOUBLE_PTR( lambda2 );

			FLA_Hev_2x2_opd( buff_alpha11,
			                 buff_alpha21,
			                 buff_alpha22,
			                 buff_lambda1,
			                 buff_lambda2 );

			break;
		}

		case FLA_COMPLEX:
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Hev_2x2_ops( float*    alpha11,
                           float*    alpha21,
                           float*    alpha22,
                           float*    lambda1,
                           float*    lambda2 )
{
	float  a11, a21, a22;
	float  l1, l2;
	float  ab, acmn, acmx, adf, df, rt, sm, tb;

	a11 = *alpha11;
	a21 = *alpha21;
	a22 = *alpha22;

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
	}
	else if ( sm > 0.0F )
	{
		l1 = 0.5F * ( sm + rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
	}
	else
	{
		l1 =  0.5F * rt;
		l2 = -0.5F * rt;
	}

	*lambda1 = l1;
	*lambda2 = l2;

	return FLA_SUCCESS;
}



FLA_Error FLA_Hev_2x2_opd( double*   alpha11,
                           double*   alpha21,
                           double*   alpha22,
                           double*   lambda1,
                           double*   lambda2 )
{
	double a11, a21, a22;
	double l1, l2;
	double ab, acmn, acmx, adf, df, rt, sm, tb;

	a11 = *alpha11;
	a21 = *alpha21;
	a22 = *alpha22;

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

	if      ( adf > ab ) rt = adf * sqrt( 1.0 + ( ab / adf ) * ( ab / adf ) );
	else if ( adf < ab ) rt = ab  * sqrt( 1.0 + ( adf / ab ) * ( adf / ab ) );
	else                 rt = ab  * sqrt( 2.0 );

	if      ( sm < 0.0 )
	{
		l1 = 0.5 * ( sm - rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
	}
	else if ( sm > 0.0 )
	{
		l1 = 0.5 * ( sm + rt );
		l2 = ( acmx / l1 ) * acmn - ( a21 / l1 ) * a21;
	}
	else
	{
		l1 =  0.5 * rt;
		l2 = -0.5 * rt;
	}

	*lambda1 = l1;
	*lambda2 = l2;

	return FLA_SUCCESS;
}

