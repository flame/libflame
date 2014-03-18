
#include "FLAME.h"

FLA_Error FLA_Gemm_hh( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_gemm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		r_val = FLA_Gemm_hh_task( alpha, A, B, beta, C, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Gemm_hh_blk_var1( alpha, A, B, beta, C, cntl );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Gemm_hh_blk_var2( alpha, A, B, beta, C, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Gemm_hh_blk_var3( alpha, A, B, beta, C, cntl );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
	{
		r_val = FLA_Gemm_hh_blk_var4( alpha, A, B, beta, C, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		r_val = FLA_Gemm_hh_blk_var5( alpha, A, B, beta, C, cntl );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT6 )
	{
		r_val = FLA_Gemm_hh_blk_var6( alpha, A, B, beta, C, cntl );
	}
#endif	
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
	{
		r_val = FLA_Gemm_hh_unb_var1( alpha, A, B, beta, C );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
	{
		r_val = FLA_Gemm_hh_unb_var2( alpha, A, B, beta, C );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
	{
		r_val = FLA_Gemm_hh_unb_var3( alpha, A, B, beta, C );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
	{
		r_val = FLA_Gemm_hh_unb_var4( alpha, A, B, beta, C );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
	{
		r_val = FLA_Gemm_hh_unb_var5( alpha, A, B, beta, C );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT6 )
	{
		r_val = FLA_Gemm_hh_unb_var6( alpha, A, B, beta, C );
	}
#endif
	else
	{
		r_val = FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

