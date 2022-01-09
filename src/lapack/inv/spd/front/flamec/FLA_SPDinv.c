/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_spdinv_t* fla_spdinv_cntl;

FLA_Error FLA_SPDinv( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_SPDinv_check( uplo, A );

  // Invoke FLA_SPDinv_internal() with an appropriate control tree.
  r_val = FLA_SPDinv_internal( uplo, A, fla_spdinv_cntl );

  return r_val;
}

