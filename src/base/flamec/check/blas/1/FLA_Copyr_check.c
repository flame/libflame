
#include "FLAME.h"

FLA_Error FLA_Copyr_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

