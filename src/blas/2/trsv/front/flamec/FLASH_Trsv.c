
#include "FLAME.h"

extern fla_trsv_t* flash_trsv_cntl;

FLA_Error FLASH_Trsv( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x )
{
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Trsv_check( uplo, trans, diag, A, x );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();

  // Execute tasks.
  r_val = FLA_Trsv_internal( uplo, trans, diag, A, x, flash_trsv_cntl );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();
  
  return r_val;
}

