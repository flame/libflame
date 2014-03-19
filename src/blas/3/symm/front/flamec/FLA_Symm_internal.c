/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_symm_t* flash_symm_cntl_blas;
extern fla_symm_t* flash_symm_cntl_mm;

FLA_Error FLA_Symm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C, fla_symm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Symm_internal_check( side, uplo, alpha, A, B, beta, C, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Symm_internal( side,
		                           uplo,
		                           alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           beta,
		                           *FLASH_OBJ_PTR_AT( C ),
		                           flash_symm_cntl_mm );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Symm( side, uplo, alpha, A, B, beta, C, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_symm_cntl_blas;
		}

		// Parameter combinations
		if      ( side == FLA_LEFT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
				r_val = FLA_Symm_ll( alpha, A, B, beta, C, cntl );
			else if ( uplo == FLA_UPPER_TRIANGULAR )
				r_val = FLA_Symm_lu( alpha, A, B, beta, C, cntl );
		}
		else if ( side == FLA_RIGHT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
				r_val = FLA_Symm_rl( alpha, A, B, beta, C, cntl );
			else if ( uplo == FLA_UPPER_TRIANGULAR )
				r_val = FLA_Symm_ru( alpha, A, B, beta, C, cntl );
		}
	}

	return r_val;
}

