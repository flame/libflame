
#include "FLAME.h"

FLA_Error FLA_Eig_gest_internal_check( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, Y );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( A, B );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
	FLA_Check_error_code( e_val );

	if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 ||
	     FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		e_val = FLA_Check_object_width_equals( Y, FLA_Obj_width( A ) );
		FLA_Check_error_code( e_val );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 ||
	          FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		e_val = FLA_Check_object_length_equals( Y, FLA_Obj_length( A ) );
		FLA_Check_error_code( e_val );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, Y, A );
		FLA_Check_error_code( e_val );
	}

	return FLA_SUCCESS;
}


