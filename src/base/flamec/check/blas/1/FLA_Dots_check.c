
#include "FLAME.h"

FLA_Error FLA_Dots_check( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj beta, FLA_Obj rho )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( x, y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( x, rho );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( x, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( x, beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( rho );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_equal_vector_dims( x, y );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

