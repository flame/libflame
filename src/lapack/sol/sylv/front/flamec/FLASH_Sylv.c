
#include "FLAME.h"

extern fla_sylv_t* flash_sylv_cntl;

FLA_Error FLASH_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale )
{
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Sylv_check( transa, transb, isgn, A, B, C, scale );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Sylv_internal( transa, transb, isgn, A, B, C, scale, flash_sylv_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  return r_val;
}

