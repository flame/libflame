/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Accum_T_UT_check( direct, storev, A, tau, T );

  // Invoke FLA_Accum_T_UT_internal().
  r_val = FLA_Accum_T_UT_internal( direct, storev, A, tau, T );

  return r_val;
}

