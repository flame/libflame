
#include "FLAME.h"

extern fla_lu_t* flash_lu_nopiv_cntl;

FLA_Error FLASH_LU_nopiv( FLA_Obj A )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_LU_nopiv_check( A );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_LU_nopiv_internal( A, flash_lu_nopiv_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  // Check for singularity.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    r_val = FLASH_LU_find_zero_on_diagonal( A );

  return r_val;
}

