/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_qrut_t*  fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT( FLA_Obj A, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_check( A, T );

  // Invoke FLA_QR_UT_internal() with the standard control tree.
  //r_val = FLA_QR_UT_internal( A, T, fla_qrut_cntl2 );
  r_val = FLA_QR_UT_internal( A, T, fla_qrut_cntl_leaf );

  return r_val;
}

