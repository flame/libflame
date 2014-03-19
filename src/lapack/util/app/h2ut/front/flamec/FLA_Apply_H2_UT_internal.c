/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_internal( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
	FLA_Error r_val = FLA_SUCCESS;

	if      ( side == FLA_LEFT )
	{
		//r_val = FLA_Apply_H2_UT_l_unb_var1( tau, u2, a1, A2 );
		r_val = FLA_Apply_H2_UT_l_opt_var1( tau, u2, a1, A2 );
	}
	else if ( side == FLA_RIGHT )
	{
		//r_val = FLA_Apply_H2_UT_r_unb_var1( tau, u2, a1, A2 );
		r_val = FLA_Apply_H2_UT_r_opt_var1( tau, u2, a1, A2 );
	}

	return r_val;
}

