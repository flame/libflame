
#include "FLAME.h"

FLA_Error FLA_CAQR_UT_inc_solve_check( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj B, FLA_Obj X )
{
  FLA_Error e_val;

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

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, X );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, ATW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, A, RTW );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, A, X, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_min( A, FLA_Obj_width( A ) * p );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

