
#include "FLAME.h"

FLA_Error FLA_LQ_UT_check( FLA_Obj A, FLA_Obj T )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( T, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

