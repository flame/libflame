
#include "FLAME.h"

FLA_Error FLA_Introduce_bulge_check( FLA_Obj shift, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon1, FLA_Obj delta2, FLA_Obj beta, FLA_Obj epsilon2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_nonconstant_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, epsilon1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( delta1, epsilon2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( gamma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( delta2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( epsilon2 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

