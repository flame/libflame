
#include "FLAME.h"

extern fla_trsv_t* flash_trsv_cntl_blas;
extern fla_trsv_t* flash_trsv_cntl;

FLA_Error FLA_Trsv_internal( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj A, FLA_Obj x, fla_trsv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Trsv_internal_check( uplo, transa, diag, A, x, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Trsv_internal( uplo,
		                           transa,
		                           diag,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( x ),
		                           flash_trsv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Trsv( uplo, transa, diag, A, x, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_trsv_cntl_blas;
		}

		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			if      ( transa == FLA_NO_TRANSPOSE )
				r_val = FLA_Trsv_ln( diag, A, x, cntl );
			else if ( transa == FLA_TRANSPOSE )
				r_val = FLA_Trsv_lt( diag, A, x, cntl );
			else if ( transa == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Trsv_lc( diag, A, x, cntl );
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			if      ( transa == FLA_NO_TRANSPOSE )
				r_val = FLA_Trsv_un( diag, A, x, cntl );
			else if ( transa == FLA_TRANSPOSE )
				r_val = FLA_Trsv_ut( diag, A, x, cntl );
			else if ( transa == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Trsv_uc( diag, A, x, cntl );
		}
	}

	return r_val;
}

