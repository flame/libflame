
#include "FLAME.h"

FLA_Error FLA_Her2k_internal_check( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_her2k_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	// Verify that the object element types are identical.
	e_val = FLA_Check_identical_object_elemtype( A, B );
	FLA_Check_error_code( e_val );

	e_val = FLA_Check_identical_object_elemtype( A, C );
	FLA_Check_error_code( e_val );

	// Verify conformality between all the objects. This check works regardless
	// of whether the element type is FLA_MATRIX or FLA_SCALAR because the
	// element length and width are used instead of scalar length and width.
	if ( trans == FLA_NO_TRANSPOSE )
	{
		e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, A, B, C );
		FLA_Check_error_code( e_val );
	}
	else
	{
		e_val = FLA_Check_matrix_matrix_dims( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, A, B, C );
		FLA_Check_error_code( e_val );
	}

	return FLA_SUCCESS;
}

