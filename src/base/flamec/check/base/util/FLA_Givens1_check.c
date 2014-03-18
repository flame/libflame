
#include "FLAME.h"

FLA_Error FLA_Givens1_check( FLA_Side side, FLA_Obj chi_1, FLA_Obj chi_2, FLA_Obj gamma, FLA_Obj sigma )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( chi_1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( chi_1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, chi_2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( chi_1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( chi_2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

