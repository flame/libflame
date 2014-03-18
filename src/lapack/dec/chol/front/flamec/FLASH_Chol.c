
#include "FLAME.h"

extern fla_chol_t* flash_chol_cntl;

FLA_Error FLASH_Chol( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Chol_check( uplo, A );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Chol_internal( uplo, A, flash_chol_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

