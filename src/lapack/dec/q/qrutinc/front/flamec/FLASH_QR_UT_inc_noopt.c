/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_qrut_t*   fla_qrut_cntl_leaf;
extern __thread fla_apqut_t*  fla_apqut_cntl_leaf;
extern __thread fla_qr2ut_t*  fla_qr2ut_cntl_leaf;
extern __thread fla_apq2ut_t* fla_apq2ut_cntl_leaf;

extern __thread fla_qrutinc_t* flash_qrutinc_cntl;

FLA_Error FLASH_QR_UT_inc_noopt( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_inc_check( A, TW );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Invoke FLA_QR_UT_inc_blk_var1() with the standard control tree.
  r_val = FLA_QR_UT_inc_blk_var1( A, TW, flash_qrutinc_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

