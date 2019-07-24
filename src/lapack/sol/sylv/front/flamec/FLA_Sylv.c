/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_sylv_t* fla_sylv_cntl;

FLA_Error FLA_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Sylv_check( transa, transb, isgn, A, B, C, scale );

  // Invoke FLA_Sylv_internal() with the appropriate control tree.
  r_val = FLA_Sylv_internal( transa, transb, isgn, A, B, C, scale, fla_sylv_cntl );

  return r_val;
}

