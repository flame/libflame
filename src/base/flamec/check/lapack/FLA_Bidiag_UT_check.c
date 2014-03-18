
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_check( FLA_Obj A, FLA_Obj TU, FLA_Obj TV )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, TU );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, TV );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( TU, FLA_Obj_length( TV ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( TU, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( TV, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

