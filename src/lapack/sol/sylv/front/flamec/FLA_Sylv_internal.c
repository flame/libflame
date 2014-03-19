/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_sylv_t* flash_sylv_cntl;
extern fla_sylv_t* fla_sylv_cntl_leaf;

FLA_Error FLA_Sylv_internal( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Sylv_internal_check( transa, transb, isgn, A, B, C, scale, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Sylv_internal( transa,
		                           transb,
		                           isgn,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           *FLASH_OBJ_PTR_AT( C ),
		                           scale,
		                           flash_sylv_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Sylv( transa, transb, isgn, A, B, C, scale, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = fla_sylv_cntl_leaf;
		}

		// Parameter combinations
		if      ( transa == FLA_NO_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Sylv_nn( isgn, A, B, C, scale, cntl );
			else if ( transb == FLA_TRANSPOSE || transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Sylv_nh( isgn, A, B, C, scale, cntl );
		}
		else if ( transa == FLA_TRANSPOSE || transa == FLA_CONJ_TRANSPOSE )
		{
			if      ( transb == FLA_NO_TRANSPOSE )
				r_val = FLA_Sylv_hn( isgn, A, B, C, scale, cntl );
			else if ( transb == FLA_TRANSPOSE || transb == FLA_CONJ_TRANSPOSE )
				r_val = FLA_Sylv_hh( isgn, A, B, C, scale, cntl );
		}
	}

	return r_val;
}

