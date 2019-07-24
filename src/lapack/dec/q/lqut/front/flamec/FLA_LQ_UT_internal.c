/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_lqut_t* flash_lqut_cntl;
extern TLS_CLASS_SPEC fla_lqut_t* flash_lqut_cntl_leaf;
extern TLS_CLASS_SPEC fla_lqut_t* fla_lqut_cntl_leaf;

FLA_Error FLA_LQ_UT_internal( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_LQ_UT_internal_check( A, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		if ( FLASH_Queue_get_enabled( ) )
		{
			// Enqueue
			ENQUEUE_FLASH_LQ_UT_macro( A, T, cntl );
		}
		else
		{
			// Execute
			r_val = FLA_LQ_UT_macro_task( A, T, cntl );
		}
	}
	else
	{
		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_LQ_UT_unb_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_LQ_UT_opt_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_LQ_UT_blk_var1( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
		{
			r_val = FLA_LQ_UT_unb_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
		{
			r_val = FLA_LQ_UT_opt_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_LQ_UT_blk_var2( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_LQ_UT_blk_var3( A, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

