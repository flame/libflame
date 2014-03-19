/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Random_symm_matrix( FLA_Uplo uplo, FLA_Obj A )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Random_symm_matrix_check( uplo, A );

  FLA_Random_tri_matrix( uplo, FLA_NONUNIT_DIAG, A );

  FLA_Symmetrize( uplo, A );

  return FLA_SUCCESS;
}

