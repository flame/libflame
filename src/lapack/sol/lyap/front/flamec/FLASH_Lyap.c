
#include "FLAME.h"

extern fla_lyap_t* flash_lyap_cntl;

FLA_Error FLASH_Lyap( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Lyap_check( trans, isgn, A, C, scale );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Lyap_internal( trans, isgn, A, C, scale, flash_lyap_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

