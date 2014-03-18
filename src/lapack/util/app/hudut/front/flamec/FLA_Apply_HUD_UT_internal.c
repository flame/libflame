
#include "FLAME.h"

FLA_Error FLA_Apply_HUD_UT_internal( FLA_Side side,
                                     FLA_Obj tau, FLA_Obj w12t,
                                                  FLA_Obj r12t,
                                     FLA_Obj u1,  FLA_Obj C2,
                                     FLA_Obj v1,  FLA_Obj D2 )
{
	FLA_Error r_val = FLA_SUCCESS;

	if      ( side == FLA_LEFT )
	{
		//r_val = FLA_Apply_HUD_UT_l_unb_var1( tau, w12t, r12t, u1, C2, v1, D2 );
		r_val = FLA_Apply_HUD_UT_l_opt_var1( tau, w12t, r12t, u1, C2, v1, D2 );
	}
	else if ( side == FLA_RIGHT )
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	return r_val;
}

