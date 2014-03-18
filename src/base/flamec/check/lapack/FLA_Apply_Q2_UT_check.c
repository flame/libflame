
#include "FLAME.h"

FLA_Error FLA_Apply_Q2_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_storev( storev );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( D, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( D, W );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( D, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( D, E );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( E );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, D, T );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_conformal_dims( FLA_TRANSPOSE, T, W );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, C, E );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, C, T );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_conformal_dims( FLA_TRANSPOSE, T, W );
    FLA_Check_error_code( e_val );

    e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, C, E );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

