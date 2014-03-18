
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_u( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
	{
		r_val = FLA_Bidiag_UT_u_unb_var1( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
	{
		r_val = FLA_Bidiag_UT_u_unb_var2( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
	{
		r_val = FLA_Bidiag_UT_u_unb_var3( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
	{
		r_val = FLA_Bidiag_UT_u_unb_var4( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
	{
		r_val = FLA_Bidiag_UT_u_unb_var5( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Bidiag_UT_u_opt_var1( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
	{
		r_val = FLA_Bidiag_UT_u_opt_var2( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
	{
		r_val = FLA_Bidiag_UT_u_opt_var3( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT4 )
	{
		r_val = FLA_Bidiag_UT_u_opt_var4( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT5 )
	{
		r_val = FLA_Bidiag_UT_u_opt_var5( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Bidiag_UT_u_blk_var1( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Bidiag_UT_u_blk_var2( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Bidiag_UT_u_blk_var3( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
	{
		r_val = FLA_Bidiag_UT_u_blk_var4( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		r_val = FLA_Bidiag_UT_u_blk_var5( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT2 )
	{
		r_val = FLA_Bidiag_UT_u_blf_var2( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT3 )
	{
		r_val = FLA_Bidiag_UT_u_blf_var3( A, TU, TV );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT4 )
	{
		r_val = FLA_Bidiag_UT_u_blf_var4( A, TU, TV );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

