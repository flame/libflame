/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Apply_pivots_ln_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	int res_cntl = FLA_Cntl_variant( cntl );

	if(res_cntl != 81)
		printf("res_cntl=%d\n", res_cntl);

	if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_EXTERN )
	{
		r_val = FLA_Apply_pivots_ln_unb_ext_ts( FLA_cntl_init_i, p, A );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Apply_pivots_ln_opt_var1_ts( FLA_cntl_init_i, p, A );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Apply_pivots_ln_blk_var1_ts( FLA_cntl_init_i, p, A, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Apply_pivots_ln_blk_var2_ts( FLA_cntl_init_i, p, A, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}
   
	return r_val;
}
#endif

FLA_Error FLA_Apply_pivots_ln( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	int res_cntl = FLA_Cntl_variant( cntl );

	if(res_cntl != 81)
		printf("res_cntl=%d\n", res_cntl);

	if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_EXTERN )
	{
		r_val = FLA_Apply_pivots_ln_unb_ext( p, A );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
	{
		r_val = FLA_Apply_pivots_ln_opt_var1( p, A );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Apply_pivots_ln_blk_var1( p, A, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Apply_pivots_ln_blk_var2( p, A, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}
   
	return r_val;
}

