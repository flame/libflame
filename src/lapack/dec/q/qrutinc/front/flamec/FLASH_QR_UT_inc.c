
#include "FLAME.h"

FLA_Error FLASH_QR_UT_inc( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;

  if ( FLASH_Queue_stack_depth() == 0 )
    r_val = FLASH_QR_UT_inc_opt1( A, TW );
  else
    r_val = FLASH_QR_UT_inc_noopt( A, TW );

  return r_val;
}

