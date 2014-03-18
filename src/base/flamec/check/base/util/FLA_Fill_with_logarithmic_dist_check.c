
#include "FLAME.h"

FLA_Error FLA_Fill_with_logarithmic_dist_check( FLA_Obj alpha, FLA_Obj x )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( x, alpha );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

