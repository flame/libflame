/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_apqut_t* flash_apqut_cntl;
extern fla_apqut_t* flash_apqut_cntl_leaf;
extern fla_apqut_t* fla_apqut_cntl_leaf;

FLA_Error FLA_Apply_Q_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Apply_Q_UT_internal_check( side, trans, direct, storev, A, T, W, B, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Apply_Q_UT_internal( side,
		                                 trans,
		                                 direct,
		                                 storev,
		                                 *FLASH_OBJ_PTR_AT( A ),
		                                 *FLASH_OBJ_PTR_AT( T ),
		                                 *FLASH_OBJ_PTR_AT( W ),
		                                 *FLASH_OBJ_PTR_AT( B ),
		                                 flash_apqut_cntl );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( A ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Apply_Q_UT( side, trans, direct, storev, A, T, W, B, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( A ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf.
			cntl = fla_apqut_cntl_leaf;
		}

		if      ( side == FLA_LEFT )
		{
			if      ( trans == FLA_NO_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_lnfc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_lnfr( A, T, W, B, cntl );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_lnbc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_lnbr( A, T, W, B, cntl );
				}
			}
			else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_lhfc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_lhfr( A, T, W, B, cntl );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_lhbc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_lhbr( A, T, W, B, cntl );
				}
			}
		}
		else if ( side == FLA_RIGHT )
		{
			if      ( trans == FLA_NO_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_rnfc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_rnfr( A, T, W, B, cntl );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_rnbc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_rnbr( A, T, W, B, cntl );
				}
			}
			else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_rhfc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_rhfr( A, T, W, B, cntl );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q_UT_rhbc( A, T, W, B, cntl );
					else if ( storev == FLA_ROWWISE )
						r_val = FLA_Apply_Q_UT_rhbr( A, T, W, B, cntl );
				}
			}
		}
	}

	return r_val;
}

