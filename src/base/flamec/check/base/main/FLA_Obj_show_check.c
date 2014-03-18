
#include "FLAME.h"

FLA_Error FLA_Obj_show_check( char* s1, FLA_Obj obj, char* format, char* s2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( s1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_scalar_elemtype( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( obj );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( format );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( s2 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

