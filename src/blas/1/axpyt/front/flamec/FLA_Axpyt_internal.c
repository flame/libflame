
#include "FLAME.h"

extern fla_axpyt_t* flash_axpyt_cntl_blas;
extern fla_axpyt_t* flash_axpyt_cntl;

FLA_Error FLA_Axpyt_internal( FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_axpyt_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Axpyt_internal_check( trans, alpha, A, B, cntl );

	if ( FLA_Obj_equals( alpha, FLA_ZERO ) )
		return FLA_SUCCESS;

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Axpyt_internal( trans,
                                    alpha,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            *FLASH_OBJ_PTR_AT( B ),
		                            flash_axpyt_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Axpyt( trans, alpha, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_axpyt_cntl_blas;
		}
		
		// Parameter combinations
		if      ( trans == FLA_NO_TRANSPOSE )
		{
			r_val = FLA_Axpyt_n( alpha, A, B, cntl );
		}
		else if ( trans == FLA_TRANSPOSE )
		{
			r_val = FLA_Axpyt_t( alpha, A, B, cntl );
		}
		else if ( trans == FLA_CONJ_NO_TRANSPOSE )
		{
			r_val = FLA_Axpyt_c( alpha, A, B, cntl );
		}
		else if ( trans == FLA_CONJ_TRANSPOSE )
		{
			r_val = FLA_Axpyt_h( alpha, A, B, cntl );
		}
	}

	return r_val;
}

