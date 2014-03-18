
#include "FLAME.h"

FLA_Error FLA_Part_2x2_check( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                          FLA_Obj *A21, FLA_Obj *A22,
                              dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A11 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A12 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A21 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A22 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_quadrant( quadrant );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

