
#include "FLAME.h"

FLA_Error FLA_Hess_UT_internal_check( FLA_Obj A, FLA_Obj T, fla_hessut_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, T );
	FLA_Check_error_code( e_val );

	return FLA_SUCCESS;
}

