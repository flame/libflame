
#include "FLAME.h"

extern fla_spdinv_t* flash_spdinv_cntl;

FLA_Error FLASH_SPDinv( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_SPDinv_check( uplo, A );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_SPDinv_internal( uplo, A, flash_spdinv_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

