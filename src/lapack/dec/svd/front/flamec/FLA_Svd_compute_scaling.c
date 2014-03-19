/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Svd_compute_scaling( FLA_Obj A, FLA_Obj sigma )
{
	FLA_Datatype dt_real;
	FLA_Obj      norm;
	FLA_Obj      safmin;
	FLA_Obj      prec;
	FLA_Obj      rmin;
	FLA_Obj      rmax;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Svd_compute_scaling_check( A, sigma );

	dt_real = FLA_Obj_datatype_proj_to_real( A );

	FLA_Obj_create( dt_real, 1, 1, 0, 0, &norm );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &prec );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &safmin );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &rmin );
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &rmax );

	// Query safmin, precision.
	FLA_Mach_params( FLA_MACH_PREC,  prec );
	FLA_Mach_params( FLA_MACH_SFMIN, safmin );

//FLA_Obj_show( "safmin", safmin, "%20.12e", "" );
//FLA_Obj_show( "prec", prec, "%20.12e", "" );

	// rmin = sqrt( safmin ) / prec;
	FLA_Copy( safmin, rmin );
	FLA_Sqrt( rmin );
	FLA_Inv_scal( prec, rmin );

	// rmax = 1 / rmin;
	FLA_Copy( rmin, rmax );
	FLA_Invert( FLA_NO_CONJUGATE, rmax );

//FLA_Obj_show( "rmin", rmin, "%20.12e", "" );
//FLA_Obj_show( "rmax", rmax, "%20.12e", "" );

	// Find the maximum absolute value of A.
	FLA_Max_abs_value( A, norm );

	if ( FLA_Obj_gt( norm, FLA_ZERO ) && FLA_Obj_lt( norm, rmin ) )
	{
		// sigma = rmin / norm;
		FLA_Copy( rmin, sigma );
		FLA_Inv_scal( norm, sigma );
	}
	else if ( FLA_Obj_gt( norm, rmax ) )
	{
		// sigma = rmax / norm;
		FLA_Copy( rmax, sigma );
		FLA_Inv_scal( norm, sigma );
	}
	else
	{
		// sigma = 1.0;
		FLA_Copy( FLA_ONE, sigma );
	}

	FLA_Obj_free( &norm );
	FLA_Obj_free( &prec );
	FLA_Obj_free( &safmin );
	FLA_Obj_free( &rmin );
	FLA_Obj_free( &rmax );

	return FLA_SUCCESS;
}

