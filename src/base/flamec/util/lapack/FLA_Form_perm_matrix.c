/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Form_perm_matrix( FLA_Obj p, FLA_Obj A )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Form_perm_matrix_check( p, A );

  // We assume that A is correctly sized, m x m, where m is the row
  // dimension of the matrix given to FLA_LU_piv() or similar function.
  FLA_Set_to_identity( A );

  // We assume that p contains pivots in native FLAME format. That is,
  // we assume the pivot type is FLA_NATIVE_PIVOTS. This is not a huge
  // assumption since the user has to go out of his way to shift the
  // pivots into LAPACK-indexed pivots.
  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );

  return FLA_SUCCESS;
}

