
#include "FLAME.h"

extern fla_copy_t* flash_copy_cntl;

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

