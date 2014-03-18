
#include "FLAME.h"

FLA_Error FLA_Ttmm_internal_check( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl )
{
	FLA_Error e_val;

	// Abort if the control structure is NULL.
	e_val = FLA_Check_null_pointer( ( void* ) cntl );
	FLA_Check_error_code( e_val );

	return FLA_SUCCESS;
}

