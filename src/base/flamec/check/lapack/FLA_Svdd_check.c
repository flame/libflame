
#include "FLAME.h"

FLA_Error FLA_Svdd_check( FLA_Svd_type jobz, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_svd_type( jobz );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( s );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, s );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, U );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, V );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( s, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( U, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( U );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( V, FLA_Obj_width( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( V );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

