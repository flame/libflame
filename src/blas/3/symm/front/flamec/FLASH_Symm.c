
#include "FLAME.h"

extern fla_symm_t* flash_symm_cntl_mm;

FLA_Error FLASH_Symm( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Symm_check( side, uplo, alpha, A, B, beta, C );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Symm_internal( side, uplo, alpha, A, B, beta, C, flash_symm_cntl_mm );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

