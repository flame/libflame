
#include "FLAME.h"

FLA_Error FLA_UDdate_UT_inc_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, C );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, D );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( R, W );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( R );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( C ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( D ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_width_equals( R, FLA_Obj_width( T ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( T, max( FLA_Obj_length( C ),
                                                  FLA_Obj_length( D ) ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( FLA_NO_TRANSPOSE, R, W );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

