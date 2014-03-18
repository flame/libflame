
#include "FLAME.h"

extern fla_apcaq2ut_t* fla_apcaq2ut_cntl_leaf;

FLA_Error FLA_Apply_CAQ2_UT_lhfc( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C,
                                                                   FLA_Obj E, fla_apcaq2ut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Apply_CAQ2_UT_lhfc_blk_var1( D, T, W, C, E, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Apply_CAQ2_UT_lhfc_blk_var2( D, T, W, C, E, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Apply_CAQ2_UT_lhfc_blk_var3( D, T, W, C, E, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

