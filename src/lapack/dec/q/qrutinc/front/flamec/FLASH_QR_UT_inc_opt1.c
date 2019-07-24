/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_qrut_t*   fla_qrut_cntl_leaf;
extern TLS_CLASS_SPEC fla_apqut_t*  fla_apqut_cntl_leaf;
extern TLS_CLASS_SPEC fla_qr2ut_t*  fla_qr2ut_cntl_leaf;
extern TLS_CLASS_SPEC fla_apq2ut_t* fla_apq2ut_cntl_leaf;

extern TLS_CLASS_SPEC fla_qrutinc_t* flash_qrutinc_cntl;

FLA_Error FLASH_QR_UT_inc_opt1( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;
  FLA_Obj   U;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_inc_check( A, TW );

  // Create a temporary matrix to hold copies of all of the blocks along the
  // diagonal of A.
  FLASH_Obj_create_diag_panel( A, &U );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Invoke FLA_QR_UT_inc_blk_var2() with the standard control tree.
  r_val = FLA_QR_UT_inc_blk_var2( A, TW, U, flash_qrutinc_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  // Free the temporary matrix.
  FLASH_Obj_free( &U );

  return r_val;
}

