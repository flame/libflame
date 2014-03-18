
#include "FLAME.h"

extern fla_copy_t* flash_copy_cntl_blas;
extern fla_copy_t* flash_copy_cntl;

FLA_Error FLA_Copy_internal( FLA_Obj A, FLA_Obj B, fla_copy_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Copy_internal_check( A, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Copy_internal( *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           flash_copy_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Copy( A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_copy_cntl_blas;
		}
		
		// Parameter combinations
		if      ( FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
		{
			r_val = FLA_Copy_task( A, B, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_Copy_blk_var1( A, B, cntl );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_Copy_blk_var2( A, B, cntl );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_Copy_blk_var3( A, B, cntl );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
		{
			r_val = FLA_Copy_blk_var4( A, B, cntl );
		}
#endif
		else
		{
			r_val = FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

