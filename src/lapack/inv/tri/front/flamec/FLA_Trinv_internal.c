/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_trinv_t* flash_trinv_cntl;
extern TLS_CLASS_SPEC fla_trinv_t* fla_trinv_cntl_leaf;

FLA_Error FLA_Trinv_internal( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Trinv_internal_check( uplo, diag, A, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Trinv_internal( uplo,
		                            diag,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            flash_trinv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Trinv( uplo, diag, A, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_trinv_cntl_leaf;
		}

		// Parameter combinations
		if      ( uplo == FLA_LOWER_TRIANGULAR )
		{
			if      ( diag == FLA_NONUNIT_DIAG )
			{
				r_val = FLA_Trinv_ln( A, cntl );
			}
			else if ( diag == FLA_UNIT_DIAG )
			{
				r_val = FLA_Trinv_lu( A, cntl );
			}
		}
		else if ( uplo == FLA_UPPER_TRIANGULAR )
		{
			if      ( diag == FLA_NONUNIT_DIAG )
			{
				r_val = FLA_Trinv_un( A, cntl );
			}
			else if ( diag == FLA_UNIT_DIAG )
			{
				r_val = FLA_Trinv_uu( A, cntl );
			}
		}
	}

	return r_val;
}

