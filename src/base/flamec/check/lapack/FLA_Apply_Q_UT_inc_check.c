
#include "FLAME.h"

FLA_Error FLA_Apply_Q_UT_inc_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B )
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

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, TW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, W1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  //e_val = FLA_Check_square( A );
  //FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, TW );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_object_width_equals( B, FLA_Obj_width( A ) );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

