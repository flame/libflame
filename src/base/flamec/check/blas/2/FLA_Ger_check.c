
#include "FLAME.h"

FLA_Error FLA_Ger_check( FLA_Obj alpha, FLA_Obj x, FLA_Obj y, FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( A, alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( y );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_vector_dims( FLA_TRANSPOSE, A, x, y );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

