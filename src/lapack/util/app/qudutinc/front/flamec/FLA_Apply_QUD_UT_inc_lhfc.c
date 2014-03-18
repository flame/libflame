
#include "FLAME.h"

FLA_Error FLA_Apply_QUD_UT_inc_lhfc( FLA_Obj T, FLA_Obj W,
                                                FLA_Obj R,
                                     FLA_Obj U, FLA_Obj C,
                                     FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl )
{
	return FLA_Apply_QUD_UT_inc_lhfc_blk_var1( T, W, R, U, C, V, D, cntl );
}

