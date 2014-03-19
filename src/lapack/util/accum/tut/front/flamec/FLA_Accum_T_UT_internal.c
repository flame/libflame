/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_internal( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T )
{
	FLA_Error r_val = FLA_SUCCESS;

	if      ( direct == FLA_FORWARD )
	{
		if      ( storev == FLA_COLUMNWISE )
			r_val = FLA_Accum_T_UT_fc_blk_var2( A, tau, T );
		else if ( storev == FLA_ROWWISE )
			r_val = FLA_Accum_T_UT_fr_blk_var2( A, tau, T );
	}
	else if ( direct == FLA_BACKWARD )
	{
		if      ( storev == FLA_COLUMNWISE )
  			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		else if ( storev == FLA_ROWWISE )
  			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

