/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Apply_G_check( side, direct, G, A );

  if ( FLA_Obj_has_zero_dim( G ) ) return FLA_SUCCESS;

  // Invoke FLA_Apply_G_internal() to parse parameters.
  r_val = FLA_Apply_G_internal( side, direct, G, A );

  return r_val;
}

