
#include "FLAME.h"

FLA_Error FLA_Apply_HUD_UT_check( FLA_Side side, FLA_Obj tau, FLA_Obj w12t, FLA_Obj r12t, FLA_Obj u1, FLA_Obj C2, FLA_Obj v1, FLA_Obj D2 )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, w12t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, r12t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, u1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, C2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, v1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( tau, D2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( w12t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( r12t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( u1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( v1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, w12t, r12t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_vector_dims( FLA_NO_TRANSPOSE, C2, r12t, u1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_vector_dims( FLA_NO_TRANSPOSE, D2, r12t, v1 );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}


