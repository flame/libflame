
#include "FLAME.h"

extern fla_copyr_t* flash_copyr_cntl;

FLA_Error FLASH_Copyr( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Copyr_check( uplo, A, B );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Execute tasks.
  r_val = FLA_Copyr_internal( uplo, A, B, flash_copyr_cntl );

  // End the parallel region.
  FLASH_Queue_end();
  
  return r_val;
}

