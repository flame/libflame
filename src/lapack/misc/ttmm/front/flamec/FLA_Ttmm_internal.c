/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_ttmm_t* flash_ttmm_cntl;
extern fla_ttmm_t* fla_ttmm_cntl_leaf;

FLA_Error FLA_Ttmm_internal( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl )
{
		FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Ttmm_internal_check( uplo, A, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Ttmm_internal( uplo,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           flash_ttmm_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Ttmm( uplo, A, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_ttmm_cntl_leaf;
		}

		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			r_val = FLA_Ttmm_l( A, cntl );
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			r_val = FLA_Ttmm_u( A, cntl );
		}
	}

	return r_val;
}

