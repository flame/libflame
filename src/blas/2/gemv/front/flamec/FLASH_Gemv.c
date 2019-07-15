/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_gemv_t* flash_gemv_cntl_fm_rp;

FLA_Error FLASH_Gemv( FLA_Trans transa, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Gemv_check( transa, alpha, A, x, beta, y );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();

  // Execute tasks.
  r_val = FLA_Gemv_internal( transa, alpha, A, x, beta, y, flash_gemv_cntl_fm_rp );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();
  
  return r_val;
}

