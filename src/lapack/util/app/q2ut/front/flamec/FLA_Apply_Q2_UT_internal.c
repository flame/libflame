/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apq2ut_t* flash_apq2ut_cntl;
extern TLS_CLASS_SPEC fla_apq2ut_t* flash_apq2ut_cntl_leaf;
extern TLS_CLASS_SPEC fla_apq2ut_t* fla_apq2ut_cntl_leaf;

FLA_Error FLA_Apply_Q2_UT_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                    FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                     FLA_Obj E, fla_apq2ut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Apply_Q2_UT_internal_check( side, trans, direct, storev, D, T, W, C, E, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( D ) == FLA_MATRIX &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		// Recurse
		r_val = FLA_Apply_Q2_UT_internal( side,
                                          trans,
                                          direct,
                                          storev,
		                                  *FLASH_OBJ_PTR_AT( D ),
		                                  *FLASH_OBJ_PTR_AT( T ),
		                                  *FLASH_OBJ_PTR_AT( W ),
		                                  *FLASH_OBJ_PTR_AT( C ),
		                                  *FLASH_OBJ_PTR_AT( E ),
		                                  flash_apq2ut_cntl_leaf );
	}
	else if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Obj_elemtype( D ) == FLA_SCALAR &&
	          FLASH_Queue_get_enabled( ) )
	{
		// Enqueue
		ENQUEUE_FLASH_Apply_Q2_UT( side, trans, direct, storev, D, T, W, C, E, cntl );
	}
	else
	{
		if ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
		     FLA_Obj_elemtype( D ) == FLA_SCALAR &&
		     !FLASH_Queue_get_enabled( ) )
		{
			// Execute leaf.
			cntl = fla_apq2ut_cntl_leaf;
		}

		if      ( side == FLA_LEFT )
		{
			if      ( trans == FLA_NO_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q2_UT_lnfc( D, T, W, C, E, cntl );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
			}
			else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
						r_val = FLA_Apply_Q2_UT_lhfc( D, T, W, C, E, cntl );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
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
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
			}
			else if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
			{
				if      ( direct == FLA_FORWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
				else if ( direct == FLA_BACKWARD )
				{
					if      ( storev == FLA_COLUMNWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
					else if ( storev == FLA_ROWWISE )
	  					FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
				}
			}
		}
	}

	return r_val;
}

