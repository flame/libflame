
#include "FLAME.h"

FLA_Error FLA_Fill_with_random_dist_check( FLA_Obj shift, FLA_Obj max, FLA_Obj x )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( shift, max );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_identical_object_precision( x, max );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( shift );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( max );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

