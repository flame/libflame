
#include "FLAME.h"

FLA_Error FLA_Tridiag_apply_Q_check( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  if ( FLA_Obj_is_real( A ) )
  {
    e_val = FLA_Check_valid_real_trans( trans );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_valid_complex_trans( trans );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, t );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_vector_dim_min( t, FLA_Obj_length( A ) - 1 );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_object_width_equals( B, FLA_Obj_length( A ) );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_vector_dim_min( t, FLA_Obj_width( A ) - 1 );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

