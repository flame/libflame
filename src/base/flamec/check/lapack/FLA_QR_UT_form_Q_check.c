
#include "FLAME.h"

FLA_Error FLA_QR_UT_form_Q_check( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, Q );
  FLA_Check_error_code( e_val );

  // The width of T is the loop guard, not A. This should not be checked.
  //e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( A ) );
  //FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( Q, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  // Q matrix should not be restricted to be a square matrix     
  // e_val = FLA_Check_square( Q );
  // FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

