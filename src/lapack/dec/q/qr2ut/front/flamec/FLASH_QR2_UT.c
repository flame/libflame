/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_qr2ut_t* flash_qr2ut_cntl;
extern fla_qr2ut_t* flash_qr2ut_cntl_leaf;
extern fla_qr2ut_t* fla_qr2ut_cntl_leaf;

FLA_Error FLASH_QR2_UT( FLA_Obj B, FLA_Obj D, FLA_Obj T )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR2_UT_check( B, D, T );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Invoke FLA_QR2_UT_internal() with the standard control tree.
  r_val = FLA_QR2_UT_internal( B, D, T, flash_qr2ut_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

