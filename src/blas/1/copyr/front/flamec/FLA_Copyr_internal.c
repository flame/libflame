/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_copyr_t* flash_copyr_cntl_blas;
extern TLS_CLASS_SPEC fla_copyr_t* flash_copyr_cntl;

FLA_Error FLA_Copyr_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, fla_copyr_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Copyr_internal_check( uplo, A, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Copyr_internal( uplo,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            *FLASH_OBJ_PTR_AT( B ),
		                            flash_copyr_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Copyr( uplo, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_copyr_cntl_blas;
		}
		
		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			r_val = FLA_Copyr_l( A, B, cntl );
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			r_val = FLA_Copyr_u( A, B, cntl );
		}
	}

	return r_val;
}

