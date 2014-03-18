
#include "FLAME.h"

FLA_Error FLA_Axpy_check( FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( A, alpha );
  FLA_Check_error_code( e_val );

  if ( FLA_Obj_is_vector( A ) && FLA_Obj_is_vector( B ) )
  {
    e_val = FLA_Check_equal_vector_dims( A, B );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, B );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

