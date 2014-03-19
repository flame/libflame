/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR2_UT_internal_check( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( B, D );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( B, T );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, B, D );
	//FLA_Check_error_code( e_val );
	e_val = FLA_Check_object_width_equals( B, FLA_Obj_width( B ) );
	FLA_Check_error_code( e_val );

	//e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, T, B, T );
	//FLA_Check_error_code( e_val );
	e_val = FLA_Check_object_width_equals( D, FLA_Obj_width( D ) );
	FLA_Check_error_code( e_val );

	return FLA_SUCCESS;
}

