
#include "FLAME.h"

FLA_Error FLA_Max_abs_value_check( FLA_Obj A, FLA_Obj amax )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

//  e_val = FLA_Check_real_object( amax );
//  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, amax );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( amax );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

