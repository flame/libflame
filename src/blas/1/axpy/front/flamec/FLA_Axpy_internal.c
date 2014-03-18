
#include "FLAME.h"

extern fla_axpy_t* flash_axpy_cntl_blas;
extern fla_axpy_t* flash_axpy_cntl;

FLA_Error FLA_Axpy_internal( FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpy_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Axpy_internal_check( alpha, A, B, cntl );

	if ( FLA_Obj_equals( alpha, FLA_ZERO ) )
		return FLA_SUCCESS;

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Axpy_internal( alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           flash_axpy_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Axpy( alpha, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_axpy_cntl_blas;
		}
		
		// Parameter combinations
		if      ( FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
		{
			r_val = FLA_Axpy_task( alpha, A, B, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_Axpy_blk_var1( alpha, A, B, cntl );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_Axpy_blk_var2( alpha, A, B, cntl );
		}
#endif
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_Axpy_blk_var3( alpha, A, B, cntl );
		}
#ifdef FLA_ENABLE_NON_CRITICAL_CODE
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT4 )
		{
			r_val = FLA_Axpy_blk_var4( alpha, A, B, cntl );
		}
#endif
		else
		{
			r_val = FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

