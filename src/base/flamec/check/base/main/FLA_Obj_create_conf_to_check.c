
#include "FLAME.h"

FLA_Error FLA_Obj_create_conf_to_check( FLA_Trans trans, FLA_Obj obj_old, FLA_Obj *obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( obj_old );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( obj );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

