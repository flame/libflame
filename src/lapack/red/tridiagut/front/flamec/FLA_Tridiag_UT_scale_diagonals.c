/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_scale_diagonals( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_scale_diagonals_check( uplo, alpha, A );
  
  if ( uplo == FLA_LOWER_TRIANGULAR )
    r_val = FLA_Bidiag_UT_l_scale_diagonals( alpha, A );
  else
    r_val = FLA_Bidiag_UT_u_scale_diagonals( alpha, A );

  return r_val;
}
