
#include "FLAME.h"

FLA_Error FLA_Obj_elem_size_check( FLA_Obj obj )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_elemtype( FLA_Obj_elemtype( obj ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_datatype( FLA_Obj_datatype( obj ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

