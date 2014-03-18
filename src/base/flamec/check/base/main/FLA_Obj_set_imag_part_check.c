
#include "FLAME.h"

FLA_Error FLA_Obj_set_imag_part_check( FLA_Obj alpha, FLA_Obj beta )
{
  FLA_Error e_val;

  e_val = FLA_Check_real_object( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( alpha, beta );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}


