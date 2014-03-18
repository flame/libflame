
#include "FLAME.h"

extern fla_trinv_t* flash_trinv_cntl;

FLA_Error FLASH_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trinv_check( uplo, diag, A );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Trinv_internal( uplo, diag, A, flash_trinv_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

