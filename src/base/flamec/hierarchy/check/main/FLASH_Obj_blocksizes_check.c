
#include "FLAME.h"

FLA_Error FLASH_Obj_blocksizes_check( FLA_Obj H, dim_t* b_m, dim_t* b_n )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( b_m );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( b_n );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

