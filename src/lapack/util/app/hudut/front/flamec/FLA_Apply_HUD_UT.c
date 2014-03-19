/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_HUD_UT( FLA_Side side,
                            FLA_Obj tau, FLA_Obj w12t,
                                         FLA_Obj r12t,
                            FLA_Obj u1,  FLA_Obj C2,
                            FLA_Obj v1,  FLA_Obj D2 )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_HUD_UT_check( side, tau, w12t, r12t, u1, C2, v1, D2 );

  // Invoke FLA_Apply_HUD_UT_internal() to parse parameters.
  r_val = FLA_Apply_HUD_UT_internal( side, tau, w12t, r12t, u1, C2, v1, D2 );

  return r_val;
}

