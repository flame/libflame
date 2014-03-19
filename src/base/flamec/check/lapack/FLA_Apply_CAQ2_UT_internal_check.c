/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_CAQ2_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apcaq2ut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( D, T );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, W );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, C );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( D, E );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects.
	if ( side == FLA_LEFT )
	{
		if ( FLA_Obj_elemtype( D ) == FLA_MATRIX )
		{
			e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, D, T );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, C, W );
			FLA_Check_error_code( e_val );
		}
		else // if ( FLA_Obj_elemtype( D ) == FLA_SCALAR )
		{
			//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, E, C, D );
			//FLA_Check_error_code( e_val );

			e_val = FLA_Check_object_width_equals( C, FLA_Obj_width( E ) );
			FLA_Check_error_code( e_val );

			e_val = FLA_Check_object_length_equals( D, FLA_Obj_length( E ) );
			FLA_Check_error_code( e_val );
		}

		//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, C, E );
		//FLA_Check_error_code( e_val );

	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return FLA_SUCCESS;
}

