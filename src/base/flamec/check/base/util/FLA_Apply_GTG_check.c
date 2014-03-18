
#include "FLAME.h"

FLA_Error FLA_Apply_GTG_check( FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon, FLA_Obj delta2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_nonconstant_object( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( gamma, sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( gamma, delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( gamma, epsilon );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( gamma, delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta2 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

