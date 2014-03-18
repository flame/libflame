
#include "FLAME.h"

FLA_Error FLA_LU_incpiv_check( FLA_Obj A, FLA_Obj p, FLA_Obj L )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, L );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_nonconstant_object( p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, p );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( L, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

