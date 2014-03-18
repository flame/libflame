
#include "FLAME.h"

FLA_Error FLA_Merge_2x1_check( FLA_Obj AT,
                               FLA_Obj AB,  FLA_Obj *A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( AT );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( AT, AB );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_adjacent_objects_2x1( AT,
                                          AB );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

