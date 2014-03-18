
#include "FLAME.h"

FLA_Error FLA_Obj_copy_view_check( FLA_Obj A, FLA_Obj* B )
{
  FLA_Error e_val;

  e_val = FLA_Check_null_pointer( A.base );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( B->base );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

