/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hess_UT_internal( FLA_Obj A, FLA_Obj T, fla_hessut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Hess_UT_internal_check( A, T, cntl );

	{
		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_Hess_UT_unb_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
		{
			r_val = FLA_Hess_UT_unb_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
		{
			r_val = FLA_Hess_UT_unb_var3( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
		{
			r_val = FLA_Hess_UT_unb_var4( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
		{
			r_val = FLA_Hess_UT_unb_var5( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_Hess_UT_opt_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
		{
			r_val = FLA_Hess_UT_opt_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
		{
			r_val = FLA_Hess_UT_opt_var3( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT4 )
		{
			r_val = FLA_Hess_UT_opt_var4( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT5 )
		{
			r_val = FLA_Hess_UT_opt_var5( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_Hess_UT_blk_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_Hess_UT_blk_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_Hess_UT_blk_var3( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
		{
			r_val = FLA_Hess_UT_blk_var4( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
		{
			r_val = FLA_Hess_UT_blk_var5( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT2 )
		{
			r_val = FLA_Hess_UT_blf_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT3 )
		{
			r_val = FLA_Hess_UT_blf_var3( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLK_FUS_VARIANT4 )
		{
			r_val = FLA_Hess_UT_blf_var4( A, T );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

