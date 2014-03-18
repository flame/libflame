
#include "FLAME.h"

FLA_Error FLA_Obj_free_buffer_check( FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( obj->base );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

