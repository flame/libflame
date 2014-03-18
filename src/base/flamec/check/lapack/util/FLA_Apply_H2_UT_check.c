
#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_check( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, u2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, a1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( u2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( a1 );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( A2, FLA_Obj_vector_dim( u2 ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_width_equals( A2, FLA_Obj_vector_dim( a1 ) );
    FLA_Check_error_code( e_val );
  }
  else // if ( side == FLA_RIGHT )
  {
    e_val = FLA_Check_object_width_equals( A2, FLA_Obj_vector_dim( u2 ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_object_length_equals( A2, FLA_Obj_vector_dim( a1 ) );
    FLA_Check_error_code( e_val );
  }

  // if columnwise 
  //e_val = FLA_Check_matrix_vector_dims( FLA_TRANSPOSE, A2, u2, a1t );
  //FLA_Check_error_code( e_val );

  // if rowwise 
  //e_val = FLA_Check_matrix_vector_dims( ... );
  //FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

