
#include "FLAME.h"

FLA_Error FLA_Obj_extract_complex_scalar_check( FLA_Obj alpha, dcomplex* alpha_value )
{
  FLA_Error e_val;

  e_val = FLA_Check_complex_object( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( alpha_value );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}


