
#include "FLAME.h"

FLA_Error FLA_QR_UT_inc_check( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, TW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, TW );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

