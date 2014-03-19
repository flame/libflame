/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_CAQ_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj R, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apcaqutinc_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( R, TW );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, W );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, B );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	if ( side == FLA_LEFT )
	{
		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, R, TW );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( R ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_width_equals( W, FLA_Obj_width( B ) );
		FLA_Check_error_code( e_val );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return FLA_SUCCESS;
}

