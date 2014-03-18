
#include "FLAME.h"

extern fla_lyap_t* flash_lyap_cntl;
extern fla_lyap_t* fla_lyap_cntl_leaf;

FLA_Error FLA_Lyap_internal( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Lyap_internal_check( trans, isgn, A, C, scale, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Lyap_internal( trans,
		                           isgn,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( C ),
		                           scale,
		                           flash_lyap_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Lyap( trans, isgn, A, C, scale, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_lyap_cntl_leaf;
		}

		// Parameter combinations
		if      ( trans == FLA_NO_TRANSPOSE )
		{
			r_val = FLA_Lyap_n( isgn, A, C, scale, cntl );
		}
		else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
		{
			r_val = FLA_Lyap_h( isgn, A, C, scale, cntl );
		}
	}

	return r_val;
}

