/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_qrut_t*  fla_qrut_piv_cntl_leaf;

FLA_Error FLA_QR_UT_piv( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_piv_check( A, T, w, p );

  FLA_Set( FLA_ZERO, w );
  FLA_QR_UT_piv_colnorm( FLA_ONE, A, w );

  r_val = FLA_QR_UT_piv_internal( A, T, w, p, fla_qrut_piv_cntl_leaf );

  return r_val;
}

