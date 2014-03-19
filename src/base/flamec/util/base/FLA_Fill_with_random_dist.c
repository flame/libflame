/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fill_with_random_dist( FLA_Obj shift, FLA_Obj max, FLA_Obj x )
{
	FLA_Obj      r, y;
	FLA_Datatype dt_real;
	dim_t        n_x;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_Fill_with_random_dist_check( shift, max, x );

	dt_real = FLA_Obj_datatype_proj_to_real( x );
	n_x     = FLA_Obj_vector_dim( x );

	FLA_Obj_create( dt_real, n_x, 1, 0, 0, &r );
	FLA_Obj_create( dt_real, n_x, 1, 0, 0, &y );

    FLA_Random_matrix( r );

    FLA_Set( FLA_ONE, y );
    FLA_Axpy( FLA_ONE, r, y );
    FLA_Scal( FLA_ONE_HALF, y );
    FLA_Scal( max, y );

    FLA_Set( shift, r );
    FLA_Axpy( FLA_ONE, y, r );

	// Overwrite x with the distribution we created in l.
	// If x is complex, then this is where the conversion between
	// datatypes happens.
	FLA_Copy( r, x );

	FLA_Obj_free( &r );
	FLA_Obj_free( &y );

	return FLA_SUCCESS;
}

