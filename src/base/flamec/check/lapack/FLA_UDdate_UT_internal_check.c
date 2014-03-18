
#include "FLAME.h"

FLA_Error FLA_UDdate_UT_internal_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( R, C );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, D );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( R, T );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_square( R );
	FLA_Check_error_code( e_val );

	if ( FLA_Obj_elemtype( R ) == FLA_MATRIX )
	{
		e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( C ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( D ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( T ) );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_object_length_equals( T, max( FLA_Obj_length( C ),
		                                                FLA_Obj_length( D ) ) );
		FLA_Check_error_code( e_val );
	}
	else
	{

	}

	return FLA_SUCCESS;
}

