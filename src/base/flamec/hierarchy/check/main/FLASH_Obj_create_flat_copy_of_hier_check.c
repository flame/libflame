
#include "FLAME.h"

FLA_Error FLASH_Obj_create_flat_copy_of_hier_check( FLA_Obj H, FLA_Obj* F )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( F );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

