
#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_check( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( T, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, w );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( w, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( p, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

