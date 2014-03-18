
#include "FLAME.h"

extern fla_scalr_t* flash_scalr_cntl;

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

