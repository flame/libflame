
#include "FLAME.h"

FLA_Error FLA_Dot_check( FLA_Obj x, FLA_Obj y, FLA_Obj rho )
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

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_equal_vector_dims( x, y );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

