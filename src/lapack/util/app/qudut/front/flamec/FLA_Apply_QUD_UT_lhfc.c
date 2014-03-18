
#include "FLAME.h"

extern fla_apqudut_t* fla_apqudut_cntl_leaf;

FLA_Error FLA_Apply_QUD_UT_lhfc( FLA_Obj T, FLA_Obj W,
                                            FLA_Obj R,
                                 FLA_Obj U, FLA_Obj C,
                                 FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if      ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
	{
		r_val = FLA_Apply_QUD_UT_lhfc_blk_var1( T, W, R, U, C, V, D, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
	{
		r_val = FLA_Apply_QUD_UT_lhfc_blk_var2( T, W, R, U, C, V, D, cntl );
	}
	else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
	{
		r_val = FLA_Apply_QUD_UT_lhfc_blk_var3( T, W, R, U, C, V, D, cntl );
	}
	else
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

