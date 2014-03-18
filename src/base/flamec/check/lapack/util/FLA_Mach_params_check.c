
#include "FLAME.h"

FLA_Error FLA_Mach_params_check( FLA_Machval machval, FLA_Obj val )
{
  FLA_Error    e_val;

  e_val = FLA_Check_valid_machval( machval );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( val );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

