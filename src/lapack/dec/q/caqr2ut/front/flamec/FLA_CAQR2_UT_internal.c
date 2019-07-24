/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_caqr2ut_t* flash_caqr2ut_cntl;
extern TLS_CLASS_SPEC fla_caqr2ut_t* fla_caqr2ut_cntl_leaf;

FLA_Error FLA_CAQR2_UT_internal( FLA_Obj B,
                                 FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_CAQR2_UT_internal_check( B, D, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( B ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_CAQR2_UT_internal( *FLASH_OBJ_PTR_AT( B ),
		                              *FLASH_OBJ_PTR_AT( D ),
		                              *FLASH_OBJ_PTR_AT( T ),
		                              flash_caqr2ut_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( B ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		if      ( FLA_Obj_structure( D ) == FLA_FULL_MATRIX )
		{
			ENQUEUE_FLASH_QR2_UT( B, D, T, cntl );
		}
		else if ( FLA_Obj_structure( D ) == FLA_UPPER_TRIANGULAR )
		{
			ENQUEUE_FLASH_CAQR2_UT( B, D, T, cntl );
		}
		else if ( FLA_Obj_structure( D ) == FLA_ZERO_MATRIX )
		{
			// Don't enqueue any tasks for zero blocks.
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( B ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf.
			if      ( FLA_Obj_structure( D ) == FLA_FULL_MATRIX )
			{
				FLA_QR2_UT_task( B, D, T, NULL );
				return FLA_SUCCESS;
			}
			else if ( FLA_Obj_structure( D ) == FLA_UPPER_TRIANGULAR )
				cntl = fla_caqr2ut_cntl_leaf;
			else if ( FLA_Obj_structure( D ) == FLA_ZERO_MATRIX )
				return FLA_SUCCESS;
			else
				FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}

		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_unb_var1( B, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_opt_var1( B, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_CAQR2_UT_blk_var1( B, D, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_CAQR2_UT_blk_var2( B, D, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

