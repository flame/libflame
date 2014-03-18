
#include "FLAME.h"

extern fla_ttmm_t* flash_ttmm_cntl;

FLA_Error FLASH_Ttmm( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Ttmm_check( uplo, A );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Ttmm_internal( uplo, A, flash_ttmm_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

