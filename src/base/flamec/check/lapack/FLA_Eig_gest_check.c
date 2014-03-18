
#include "FLAME.h"

FLA_Error FLA_Eig_gest_check( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_inverse( inv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( B );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

