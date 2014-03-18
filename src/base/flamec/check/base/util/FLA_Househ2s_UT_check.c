
#include "FLAME.h"

FLA_Error FLA_Househ2s_UT_check( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj alpha, FLA_Obj chi_1_minus_alpha, FLA_Obj tau )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( chi_1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, x2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, chi_1_minus_alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( chi_1, tau );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( chi_1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( chi_1_minus_alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( tau );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

