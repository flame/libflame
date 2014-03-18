
#include "FLAME.h"

FLA_Error FLASH_Obj_create_helper_check( FLA_Bool without_buffer, FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( b_m );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( b_n );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( H );
  FLA_Check_error_code( e_val );

  // A value of depth < 0 should cause an error.

  // Values of m < 1, n < 1 should cause an error. (or < 0?)

  // First depth entries in blocksize_m, _n should be checked; values < 1 should cause error.

  return FLA_SUCCESS;
}

