/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_chol_t* fla_chol_cntl;
extern __thread fla_chol_t* fla_chol_cntl2;

FLA_Error FLA_Chol( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Chol_check( uplo, A );

  // Invoke FLA_Chol_internal() with the appropriate control tree.
  r_val = FLA_Chol_internal( uplo, A, fla_chol_cntl2 );
  //r_val = FLA_Chol_internal( uplo, A, fla_chol_cntl );

  return r_val;
}

