
#include "FLAME.h"

FLA_Error FLA_Obj_extract_real_scalar_check( FLA_Obj alpha, double* alpha_value )
{
  FLA_Error e_val;

  e_val = FLA_Check_real_object( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( alpha_value );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}


