/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_trsm_t* flash_trsm_cntl_blas;
extern __thread fla_trsm_t* flash_trsm_cntl_mm;

FLA_Error FLA_Trsm_internal( FLA_Side side, FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, fla_trsm_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Trsm_internal_check( side, uplo, transa, diag, alpha, A, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Trsm_internal( side,
		                           uplo,
		                           transa,
		                           diag,
		                           alpha,
		                           *FLASH_OBJ_PTR_AT( A ),
		                           *FLASH_OBJ_PTR_AT( B ),
		                           flash_trsm_cntl_mm );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Trsm( side, uplo, transa, diag, alpha, A, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf
			cntl = flash_trsm_cntl_blas;
		}

		// Parameter combinations
		if      ( side == FLA_LEFT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trsm_lln( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trsm_llt( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trsm_llc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trsm_llh( diag, alpha, A, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trsm_lun( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trsm_lut( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trsm_luc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trsm_luh( diag, alpha, A, B, cntl );
			}
		}
		else if ( side == FLA_RIGHT )
		{
			if      ( uplo == FLA_LOWER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trsm_rln( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trsm_rlt( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trsm_rlc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trsm_rlh( diag, alpha, A, B, cntl );
			}
			else if ( uplo == FLA_UPPER_TRIANGULAR )
			{
				if      ( transa == FLA_NO_TRANSPOSE )
					r_val = FLA_Trsm_run( diag, alpha, A, B, cntl );
				else if ( transa == FLA_TRANSPOSE )
					r_val = FLA_Trsm_rut( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_NO_TRANSPOSE )
					r_val = FLA_Trsm_ruc( diag, alpha, A, B, cntl );
				else if ( transa == FLA_CONJ_TRANSPOSE )
					r_val = FLA_Trsm_ruh( diag, alpha, A, B, cntl );
			}
		}
	}

	return r_val;
}

