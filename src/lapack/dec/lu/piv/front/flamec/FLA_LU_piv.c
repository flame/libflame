/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_lu_t* fla_lu_piv_cntl;
extern fla_lu_t* fla_lu_piv_cntl2;

FLA_Error FLA_LU_piv( FLA_Obj A, FLA_Obj p )
{
  FLA_Error r_val = FLA_SUCCESS;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  // Invoke FLA_LU_piv_internal() with large control tree.
  r_val = FLA_LU_piv_internal( A, p, fla_lu_piv_cntl2 );

  // This is invalid as FLA_LU_piv_internal returns a null pivot index.
  // Check for singularity.
  //if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  //  r_val = FLA_LU_find_zero_on_diagonal( A );

  return r_val;
}

