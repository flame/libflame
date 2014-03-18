
#include "FLAME.h"

FLA_Error FLA_Apply_CAQ_UT_inc_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj W, FLA_Obj B )
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

  e_val = FLA_Check_identical_object_datatype( A, ATW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, RTW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, W );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, ATW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, RTW );
  FLA_Check_error_code( e_val );

  if ( side == FLA_LEFT )
  {
    e_val = FLA_Check_object_length_equals( B, FLA_Obj_length( A ) );
    FLA_Check_error_code( e_val );
  }
  else
  {
    FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
  }

  return FLA_SUCCESS;
}

