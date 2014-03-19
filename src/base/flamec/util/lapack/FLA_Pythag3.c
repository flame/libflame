/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Pythag3( FLA_Obj chi, FLA_Obj psi, FLA_Obj zeta, FLA_Obj rho )
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( chi );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  buff_chi  = FLA_FLOAT_PTR( chi );
			float*  buff_psi  = FLA_FLOAT_PTR( psi );
			float*  buff_zeta = FLA_FLOAT_PTR( zeta );
			float*  buff_rho  = FLA_FLOAT_PTR( rho );

			FLA_Pythag3_ops( buff_chi,
			                 buff_psi,
			                 buff_zeta,
			                 buff_rho );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_chi  = FLA_DOUBLE_PTR( chi );
			double* buff_psi  = FLA_DOUBLE_PTR( psi );
			double* buff_zeta = FLA_DOUBLE_PTR( zeta );
			double* buff_rho  = FLA_DOUBLE_PTR( rho );

			FLA_Pythag3_opd( buff_chi,
			                 buff_psi,
			                 buff_zeta,
			                 buff_rho );

			break;
		}

		case FLA_COMPLEX:
		{
			FLA_Check_error_code( FLA_OBJECT_NOT_REAL );

			break;
		}

		case FLA_DOUBLE_COMPLEX:
		{
			FLA_Check_error_code( FLA_OBJECT_NOT_REAL );

			break;
		}
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Pythag3_ops( float*    chi,
                           float*    psi,
                           float*    zeta,
                           float*    rho )
{
	float  zero = bl1_s0();

	float  xabs, yabs, zabs;
	float  w;
	float  xabsdivw;
	float  yabsdivw;
	float  zabsdivw;

	xabs = fabsf( *chi );
	yabs = fabsf( *psi );
	zabs = fabsf( *zeta );
	w    = max( xabs, max( yabs, zabs ) );

	if ( w == zero )
	{
		// From netlib dlapy3:
		// W can be zero for max(0,nan,0). Adding all three entries
		// together will make sure NaN will not disappear.
		*rho = xabs + yabs + zabs;
	}
	else
	{
		xabsdivw = xabs / w;
		yabsdivw = yabs / w;
		zabsdivw = zabs / w;

		*rho = w * sqrt( xabsdivw * xabsdivw +
		                 yabsdivw * yabsdivw +
		                 zabsdivw * zabsdivw );
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Pythag3_opd( double*   chi,
                           double*   psi,
                           double*   zeta,
                           double*   rho )
{
	double zero = bl1_d0();

	double xabs, yabs, zabs;
	double w;
	double xabsdivw;
	double yabsdivw;
	double zabsdivw;

	xabs = fabs( *chi );
	yabs = fabs( *psi );
	zabs = fabs( *zeta );
	w    = max( xabs, max( yabs, zabs ) );

	if ( w == zero )
	{
		// From netlib dlapy3:
		// W can be zero for max(0,nan,0). Adding all three entries
		// together will make sure NaN will not disappear.
		*rho = xabs + yabs + zabs;
	}
	else
	{
		xabsdivw = xabs / w;
		yabsdivw = yabs / w;
		zabsdivw = zabs / w;

		*rho = w * sqrt( xabsdivw * xabsdivw +
		                 yabsdivw * yabsdivw +
		                 zabsdivw * zabsdivw );
	}

	return FLA_SUCCESS;
}

