
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_recover_tau_check( FLA_Obj TU, FLA_Obj TV, FLA_Obj tu, FLA_Obj tv )
{
  FLA_Error    e_val;

  e_val = FLA_Check_floating_object( TU );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( TU, TV );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( TU, tu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( TU, tv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( tu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( tv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( TU, FLA_Obj_width( TV ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( TU, FLA_Obj_vector_dim( tu ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( TV, FLA_Obj_vector_dim( tv ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

