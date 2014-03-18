
#include "FLAME.h"

FLA_Error FLA_Lyap_internal_check( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, C );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
    e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, C );
    FLA_Check_error_code( e_val );

	e_val = FLA_Check_square( A );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_square( C );
	FLA_Check_error_code( e_val );

	return FLA_SUCCESS;
}

