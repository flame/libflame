
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_l( FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
	{
		r_val = FLA_Tridiag_UT_l_unb_var1( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
	{
		r_val = FLA_Tridiag_UT_l_unb_var2( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
	{
		r_val = FLA_Tridiag_UT_l_unb_var3( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Tridiag_UT_l_opt_var1( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
	{
		r_val = FLA_Tridiag_UT_l_opt_var2( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
	{
		r_val = FLA_Tridiag_UT_l_opt_var3( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Tridiag_UT_l_blk_var1( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Tridiag_UT_l_blk_var2( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Tridiag_UT_l_blk_var3( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT2 )
	{
		r_val = FLA_Tridiag_UT_l_blf_var2( A, T );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT3 )
	{
		r_val = FLA_Tridiag_UT_l_blf_var3( A, T );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

