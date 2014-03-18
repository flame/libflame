
#include "FLAME.h"

FLA_Error FLA_Trsv_ln( FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		r_val = FLA_Trsv_ln_task( diag, A, x, cntl );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Trsv_ln_blk_var1( diag, A, x, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Trsv_ln_blk_var2( diag, A, x, cntl );
	}
	else
	{
		r_val = FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

