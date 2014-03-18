
#include "FLAME.h"

extern fla_uddateut_t* flash_uddateut_cntl;
extern fla_uddateut_t* fla_uddateut_cntl_leaf;

FLA_Error FLA_UDdate_UT_internal( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_UDdate_UT_internal_check( R, C, D, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( R ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_UDdate_UT_internal( *FLASH_OBJ_PTR_AT( R ),
		                                *FLASH_OBJ_PTR_AT( C ),
		                                *FLASH_OBJ_PTR_AT( D ),
		                                *FLASH_OBJ_PTR_AT( T ),
		                                flash_uddateut_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( R ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_UDdate_UT( R, C, D, T, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( R ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf.
			cntl = fla_uddateut_cntl_leaf;
		}

		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_UDdate_UT_unb_var1( R, C, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_UDdate_UT_opt_var1( R, C, D, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_UDdate_UT_blk_var1( R, C, D, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_UDdate_UT_blk_var2( R, C, D, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

