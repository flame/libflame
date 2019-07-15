/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_lu_t* fla_lu_nopiv_cntl;
extern __thread fla_lu_t* fla_lu_nopiv_cntl2;

FLA_Error FLA_LU_nopiv( FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_nopiv_check( A );

  // Invoke FLA_LU_nopiv_internal() with large control tree.
  r_val = FLA_LU_nopiv_internal( A, fla_lu_nopiv_cntl2 );

  // Check for singularity.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    r_val = FLA_LU_find_zero_on_diagonal( A );

  return r_val;
}

