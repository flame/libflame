/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_gemv_t* flash_gemv_cntl_blas;
extern TLS_CLASS_SPEC fla_gemv_t* flash_gemv_cntl_fm_rp;

FLA_Error FLA_Gemv_internal( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y, fla_gemv_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Gemv_internal_check( transa, alpha, A, x, beta, y, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Gemv_internal( transa, 
		                           alpha, 
		                           *FLASH_OBJ_PTR_AT( A ), 
		                           *FLASH_OBJ_PTR_AT( x ), 
		                           beta, 
		                           *FLASH_OBJ_PTR_AT( y ), 
		                           flash_gemv_cntl_fm_rp );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Gemv( transa, alpha, A, x, beta, y, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_gemv_cntl_blas;
		}

		// Parameter combinations
		if      ( transa == FLA_NO_TRANSPOSE )
		{
			r_val = FLA_Gemv_n( alpha, A, x, beta, y, cntl );
		}
		else if ( transa == FLA_TRANSPOSE )
		{
			r_val = FLA_Gemv_t( alpha, A, x, beta, y, cntl );
		}
		else if ( transa == FLA_CONJ_TRANSPOSE )
		{
			r_val = FLA_Gemv_h( alpha, A, x, beta, y, cntl );
		}
	}

	return r_val;
}

