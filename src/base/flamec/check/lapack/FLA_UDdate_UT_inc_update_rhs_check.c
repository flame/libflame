
#include "FLAME.h"

FLA_Error FLA_UDdate_UT_inc_update_rhs_check( FLA_Obj T, FLA_Obj bR, FLA_Obj C, FLA_Obj bC, FLA_Obj D, FLA_Obj bD )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bR );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bC );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( T, bD );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( T, max( FLA_Obj_length( C ),
                                                  FLA_Obj_length( D ) ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( C, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( D, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, C, bR, bC );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_matrix_matrix_dims( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, D, bR, bD );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

