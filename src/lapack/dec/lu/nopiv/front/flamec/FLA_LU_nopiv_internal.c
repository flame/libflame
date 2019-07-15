/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_lu_t* flash_lu_nopiv_cntl;
extern __thread fla_lu_t* fla_lu_nopiv_cntl_leaf;

FLA_Error FLA_LU_nopiv_internal( FLA_Obj A, fla_lu_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_LU_nopiv_internal_check( A, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_LU_nopiv_internal( *FLASH_OBJ_PTR_AT( A ),
		                               flash_lu_nopiv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_LU_nopiv( A, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_lu_nopiv_cntl_leaf;
		}
		
		// Choose implementation.
		if      ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_LU_nopiv_opt_var1( A );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
		{
			r_val = FLA_LU_nopiv_opt_var2( A );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT3 )
		{
			r_val = FLA_LU_nopiv_opt_var3( A );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT4 )
		{
			r_val = FLA_LU_nopiv_opt_var4( A );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT5 )
		{
			r_val = FLA_LU_nopiv_opt_var5( A );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_LU_nopiv_unb_var1( A );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
		{
			r_val = FLA_LU_nopiv_unb_var2( A );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT3 )
		{
			r_val = FLA_LU_nopiv_unb_var3( A );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT4 )
		{
			r_val = FLA_LU_nopiv_unb_var4( A );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT5 )
		{
			r_val = FLA_LU_nopiv_unb_var5( A );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_LU_nopiv_blk_var1( A, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_LU_nopiv_blk_var2( A, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_LU_nopiv_blk_var3( A, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
		{
			r_val = FLA_LU_nopiv_blk_var4( A, cntl );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT5 )
		{
			r_val = FLA_LU_nopiv_blk_var5( A, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

