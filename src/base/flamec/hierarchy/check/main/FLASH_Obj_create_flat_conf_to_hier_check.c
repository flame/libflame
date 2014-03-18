
#include "FLAME.h"

FLA_Error FLASH_Obj_create_flat_conf_to_hier_check( FLA_Trans trans, FLA_Obj H, FLA_Obj* F )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_real_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( F );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

