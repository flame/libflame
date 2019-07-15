/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_scalr_t* flash_scalr_cntl_blas;
extern __thread fla_scalr_t* flash_scalr_cntl;

FLA_Error FLA_Scalr_internal( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, fla_scalr_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Scalr_internal_check( uplo, alpha, A, cntl );

	if ( FLA_Obj_equals( alpha, FLA_ONE ) )
		return FLA_SUCCESS;

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Scalr_internal( uplo,
		                            alpha,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            flash_scalr_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Scalr( uplo, alpha, A, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_scalr_cntl_blas;
		}
		
		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			r_val = FLA_Scalr_l( alpha, A, cntl );
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			r_val = FLA_Scalr_u( alpha, A, cntl );
		}
	}

	return r_val;
}

