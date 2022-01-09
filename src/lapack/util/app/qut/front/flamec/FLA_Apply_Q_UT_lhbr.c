/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_apqut_t* fla_apqut_cntl_leaf;

FLA_Error FLA_Apply_Q_UT_lhbr( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Apply_Q_UT_lhbr_blk_var1( A, T, W, B, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Apply_Q_UT_lhbr_blk_var2( A, T, W, B, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Apply_Q_UT_lhbr_blk_var3( A, T, W, B, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

