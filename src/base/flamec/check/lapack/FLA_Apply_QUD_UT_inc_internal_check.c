/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( R, T );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, W );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, U );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, C );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, V );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, D );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects.
	if ( side == FLA_LEFT )
	{
		e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( U ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_length_equals( T, fla_max( FLA_Obj_length( U ),
		                                                FLA_Obj_length( V ) ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, W, R );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, U, R, C );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, V, R, D );
		FLA_Check_error_code( e_val );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return FLA_SUCCESS;
}

