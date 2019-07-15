/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_copy_t* flash_copy_cntl;

FLA_Error FLASH_Copy( FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Copy_check( A, B );

  // Begin a parallel region.
  FLASH_Queue_begin();

  // Execute tasks.
  r_val = FLA_Copy_internal( A, B, flash_copy_cntl );

  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

