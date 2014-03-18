
#include "FLAME.h"

FLA_Error FLA_Max_abs_value_herm_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj maxabs )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( maxabs );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( maxabs );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, maxabs );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( maxabs );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

