/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_nn( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_EXTERN )
	{
		r_val = FLA_Sylv_nn_blk_ext( isgn, A, B, C, scale );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_EXTERN )
	{
		r_val = FLA_Sylv_nn_unb_ext( isgn, A, B, C, scale );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Sylv_nn_opt_var1( isgn, A, B, C, scale );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Sylv_nn_blk_var1( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Sylv_nn_blk_var2( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Sylv_nn_blk_var3( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
	{
		r_val = FLA_Sylv_nn_blk_var4( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
	{
		r_val = FLA_Sylv_nn_blk_var5( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT6 )
	{
		r_val = FLA_Sylv_nn_blk_var6( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT7 )
	{
		r_val = FLA_Sylv_nn_blk_var7( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT8 )
	{
		r_val = FLA_Sylv_nn_blk_var8( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT9 )
	{
		r_val = FLA_Sylv_nn_blk_var9( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT10 )
	{
		r_val = FLA_Sylv_nn_blk_var10( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT11 )
	{
		r_val = FLA_Sylv_nn_blk_var11( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT12 )
	{
		r_val = FLA_Sylv_nn_blk_var12( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT13 )
	{
		r_val = FLA_Sylv_nn_blk_var13( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT14 )
	{
		r_val = FLA_Sylv_nn_blk_var14( isgn, A, B, C, scale, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT15 )
	{
		r_val = FLA_Sylv_nn_blk_var15( isgn, A, B, C, scale, cntl );
	}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT16 )
	{
		r_val = FLA_Sylv_nn_blk_var16( isgn, A, B, C, scale, cntl );
	}
#endif
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT17 )
	{
		r_val = FLA_Sylv_nn_blk_var17( isgn, A, B, C, scale, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT18 )
	{
		r_val = FLA_Sylv_nn_blk_var18( isgn, A, B, C, scale, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

