/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_appiv_t* fla_appiv_cntl_leaf;

FLA_Error FLA_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  
  r_val = FLA_Apply_pivots_internal( side, trans, p, A, fla_appiv_cntl_leaf );
  
  return r_val;
}

