
#include "FLAME.h"

FLA_Error FLA_Eig_gest_nl( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_EXTERN )
	{
		r_val = FLA_Eig_gest_nl_blk_ext( A, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_EXTERN )
	{
		r_val = FLA_Eig_gest_nl_unb_ext( A, B );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
	{
		r_val = FLA_Eig_gest_nl_unb_var1( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
	{
		r_val = FLA_Eig_gest_nl_unb_var2( A, Y, B );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
	{
		//r_val = FLA_Eig_gest_nl_unb_var3( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
	{
		r_val = FLA_Eig_gest_nl_unb_var4( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
	{
		r_val = FLA_Eig_gest_nl_unb_var5( A, Y, B );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Eig_gest_nl_opt_var1( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
	{
		r_val = FLA_Eig_gest_nl_opt_var2( A, Y, B );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
	{
		//r_val = FLA_Eig_gest_nl_opt_var3( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT4 )
	{
		r_val = FLA_Eig_gest_nl_opt_var4( A, Y, B );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT5 )
	{
		r_val = FLA_Eig_gest_nl_opt_var5( A, Y, B );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Eig_gest_nl_blk_var1( A, Y, B, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Eig_gest_nl_blk_var2( A, Y, B, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		//r_val = FLA_Eig_gest_nl_blk_var3( A, Y, B, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
	{
		r_val = FLA_Eig_gest_nl_blk_var4( A, Y, B, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		r_val = FLA_Eig_gest_nl_blk_var5( A, Y, B, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

