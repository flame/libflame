
#include "FLAME.h"

FLA_Error FLA_Merge_2x2_check( FLA_Obj A11, FLA_Obj A12,
                               FLA_Obj A21, FLA_Obj A22,  FLA_Obj *A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( A11, A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( A12, A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_base_buffer_mismatch( A11, A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_adjacent_objects_2x2( A11, A12,
                                          A21, A22 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

