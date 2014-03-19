/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Pythag2( FLA_Obj chi, FLA_Obj psi, FLA_Obj rho )
{
	FLA_Datatype datatype;

	datatype = FLA_Obj_datatype( chi );

	switch ( datatype )
	{
		case FLA_FLOAT:
		{
			float*  buff_chi = FLA_FLOAT_PTR( chi );
			float*  buff_psi = FLA_FLOAT_PTR( psi );
			float*  buff_rho = FLA_FLOAT_PTR( rho );

			FLA_Pythag2_ops( buff_chi,
			                 buff_psi,
			                 buff_rho );

			break;
		}

		case FLA_DOUBLE:
		{
			double* buff_chi = FLA_DOUBLE_PTR( chi );
			double* buff_psi = FLA_DOUBLE_PTR( psi );
			double* buff_rho = FLA_DOUBLE_PTR( rho );

			FLA_Pythag2_opd( buff_chi,
			                 buff_psi,
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



FLA_Error FLA_Pythag2_ops( float*    chi,
                           float*    psi,
                           float*    rho )
{
	float  zero = bl1_s0();
	float  one  = bl1_s1();

	float  xabs, yabs;
	float  w, z;
	float  zdivw;

	xabs = fabsf( *chi );
	yabs = fabsf( *psi );
	w    = max( xabs, yabs );
	z    = min( xabs, yabs );

	if ( z == zero )
	{
		*rho = w;
	}
	else
	{
		zdivw = z / w;

		*rho = w * sqrt( one + zdivw * zdivw );
	}

	return FLA_SUCCESS;
}



FLA_Error FLA_Pythag2_opd( double*   chi,
                           double*   psi,
                           double*   rho )
{
	double zero = bl1_d0();
	double one  = bl1_d1();

	double xabs, yabs;
	double w, z;
	double zdivw;

	xabs = fabs( *chi );
	yabs = fabs( *psi );
	w    = max( xabs, yabs );
	z    = min( xabs, yabs );

	if ( z == zero )
	{
		*rho = w;
	}
	else
	{
		zdivw = z / w;

		*rho = w * sqrt( one + zdivw * zdivw );
	}

	return FLA_SUCCESS;
}

