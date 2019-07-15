/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_tridiagut_t* fla_tridiagut_cntl_fused;
extern __thread fla_tridiagut_t* fla_tridiagut_cntl_nofus;
extern __thread fla_tridiagut_t* fla_tridiagut_cntl_plain;

FLA_Error FLA_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_check( uplo, A, T );

  if ( FLA_Obj_row_stride( A ) == 1 &&
       FLA_Obj_is_double_precision( A ) )
    // Temporary fix not to use the fused version (numerically unstable).
    r_val = FLA_Tridiag_UT_internal( uplo, A, T, fla_tridiagut_cntl_plain );
  else
    r_val = FLA_Tridiag_UT_internal( uplo, A, T, fla_tridiagut_cntl_nofus );

  return r_val;
}

