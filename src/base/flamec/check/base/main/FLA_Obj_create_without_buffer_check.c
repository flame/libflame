
#include "FLAME.h"

FLA_Error FLA_Obj_create_without_buffer_check( FLA_Datatype datatype, dim_t m, dim_t n, FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_datatype( datatype );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

