/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_copyt_t* flash_copyt_cntl_blas;
extern TLS_CLASS_SPEC fla_copyt_t* flash_copyt_cntl;

FLA_Error FLA_Copyt_internal( FLA_Trans trans, FLA_Obj A, FLA_Obj B, fla_copyt_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Copyt_internal_check( trans, A, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Copyt_internal( trans,
		                            *FLASH_OBJ_PTR_AT( A ),
		                            *FLASH_OBJ_PTR_AT( B ),
		                            flash_copyt_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Copyt( trans, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_copyt_cntl_blas;
		}
		
		// Parameter combinations
		if      ( trans == FLA_NO_TRANSPOSE )
		{
			r_val = FLA_Copyt_n( A, B, cntl );
		}
		else if ( trans == FLA_TRANSPOSE )
		{
			r_val = FLA_Copyt_t( A, B, cntl );
		}
		else if ( trans == FLA_CONJ_NO_TRANSPOSE )
		{
			r_val = FLA_Copyt_c( A, B, cntl );
		}
		else if ( trans == FLA_CONJ_TRANSPOSE )
		{
			r_val = FLA_Copyt_h( A, B, cntl );
		}
	}

	return r_val;
}

