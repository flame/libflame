
#include "FLAME.h"

extern fla_trmm_t* flash_trmm_cntl_mm;

FLA_Error FLASH_Trmm( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trmm_check( side, uplo, trans, diag, alpha, A, B );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Trmm_internal( side, uplo, trans, diag, alpha, A, B, flash_trmm_cntl_mm );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

