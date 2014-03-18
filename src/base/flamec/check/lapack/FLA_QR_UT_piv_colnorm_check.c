
#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_colnorm_check(  FLA_Obj alpha, FLA_Obj A, FLA_Obj b )
{
  FLA_Error e_val;

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, b );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( b, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

