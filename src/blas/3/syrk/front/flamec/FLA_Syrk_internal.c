/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_syrk_t* flash_syrk_cntl_blas;
extern __thread fla_syrk_t* flash_syrk_cntl_mm;

FLA_Error FLA_Syrk_internal( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C, fla_syrk_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Syrk_internal_check( uplo, trans, alpha, A, beta, C, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Syrk_internal( uplo,
		                           trans,
		                           alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           beta,
		                           *FLASH_OBJ_PTR_AT( C ),
		                           flash_syrk_cntl_mm );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Syrk( uplo, trans, alpha, A, beta, C, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_syrk_cntl_blas;
		}

		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			if      ( trans == FLA_NO_TRANSPOSE )
				r_val = FLA_Syrk_ln( alpha, A, beta, C, cntl );
			else if ( trans == FLA_TRANSPOSE )
				r_val = FLA_Syrk_lt( alpha, A, beta, C, cntl );
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			if      ( trans == FLA_NO_TRANSPOSE )
				r_val = FLA_Syrk_un( alpha, A, beta, C, cntl );
			else if ( trans == FLA_TRANSPOSE )
				r_val = FLA_Syrk_ut( alpha, A, beta, C, cntl );
		}
	}

	return r_val;
}

