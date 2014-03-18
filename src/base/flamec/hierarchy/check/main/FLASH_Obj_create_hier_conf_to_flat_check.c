
#include "FLAME.h"

FLA_Error FLASH_Obj_create_hier_conf_to_flat_check( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* blocksizes, FLA_Obj* H )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_real_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( blocksizes );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( H );
  FLA_Check_error_code( e_val );

  // A value of depth < 0 should cause an error.

  // First depth entries in blocksize should be checked; values < 1 should cause error.

  return FLA_SUCCESS;
}

