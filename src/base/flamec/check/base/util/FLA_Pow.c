
#include "FLAME.h"

FLA_Error FLA_Pow_check( FLA_Obj base, FLA_Obj exp, FLA_Obj btoe )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( base );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( exp );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( btoe );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( base );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( exp );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( btoe );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( btoe );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

