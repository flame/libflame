
#include "FLAME.h"

FLA_Error FLA_Part_1x2_check( FLA_Obj A,  FLA_Obj *A1, FLA_Obj *A2,
                                          dim_t    nb, FLA_Side side )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A1 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( A2 );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_leftright_side( side );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

