
#include "FLAME.h"

FLA_Error FLA_SPDinv_internal( FLA_Uplo uplo, FLA_Obj A, fla_spdinv_t* cntl )
{
	FLA_Error r_val;
	FLA_Error e_val;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_SPDinv_internal_check( uplo, A, cntl );

	r_val = FLA_Chol_internal( uplo, A,
	                           FLA_Cntl_sub_chol( cntl ) );

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_chol_failure( r_val );
		FLA_Check_error_code( e_val );
	}

	FLA_Trinv_internal( uplo, FLA_NONUNIT_DIAG, A, 
	                    FLA_Cntl_sub_trinv( cntl ) );

	FLA_Ttmm_internal( uplo, A, 
	                   FLA_Cntl_sub_ttmm( cntl ) );

	return FLA_SUCCESS;
}

