
#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, fla_apqutinc_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, TW );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( A, W1 );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( A, B );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	if ( side == FLA_LEFT )
	{
		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, TW );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_width_equals( W1, FLA_Obj_width( B ) );
		FLA_Check_error_code( e_val );
	}
	else
	{
		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, TW );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_width_equals( B, FLA_Obj_width( A ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_length_equals( W1, FLA_Obj_length( B ) );
		FLA_Check_error_code( e_val );
	}

	return FLA_SUCCESS;
}

