/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_CAQ_UT_inc_internal( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev,
                                         FLA_Obj R, FLA_Obj TW, FLA_Obj W, FLA_Obj B,
                                         fla_apcaqutinc_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Apply_CAQ_UT_inc_internal_check( side, trans, direct, storev, R, TW, W, B, cntl );

	if      ( side == FLA_LEFT )
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
					r_val = FLA_Apply_CAQ_UT_inc_lhfc( R, TW, W, B, cntl );
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

	return r_val;
}

