
#include "FLAME.h"

FLA_Error FLA_Hevd_compute_scaling_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj sigma )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( sigma );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, sigma );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( sigma );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

