/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_scalr_t* flash_scalr_cntl;

FLA_Error FLASH_Scalr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Scalr_check( uplo, alpha, A );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();

  // Execute tasks.
  r_val = FLA_Scalr_internal( uplo, alpha, A, flash_scalr_cntl );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();
  
  return r_val;
}

