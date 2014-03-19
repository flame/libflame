/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_tpose_t* fla_tpose_cntl;

FLA_Error FLA_Transpose( FLA_Obj A )
{
  FLA_Error r_val;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Transpose_check( A );

  r_val = FLA_Transpose_blk_var2( A, fla_tpose_cntl );

  return r_val;
}

