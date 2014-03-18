
#include "FLAME.h"

extern fla_axpy_t* flash_axpy_cntl;

FLA_Error FLASH_Axpy( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error r_val;
  FLA_Bool  enable_supermatrix;
  
  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Axpy_check( alpha, A, B );

  // Find the status of SuperMatrix.
  enable_supermatrix = FLASH_Queue_get_enabled();

  // Temporarily disable SuperMatrix.
  FLASH_Queue_disable();

  // Execute tasks.
  r_val = FLA_Axpy_internal( alpha, A, B, flash_axpy_cntl );

  // Restore SuperMatrix to its previous status.
  if ( enable_supermatrix )
     FLASH_Queue_enable();
  
  return r_val;
}

