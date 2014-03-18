
#include "FLAME.h"

extern fla_her2k_t* flash_her2k_cntl_mm;

FLA_Error FLASH_Her2k( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Error r_val;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Her2k_check( uplo, trans, alpha, A, B, beta, C );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Her2k_internal( uplo, trans, alpha, A, B, beta, C, flash_her2k_cntl_mm );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

