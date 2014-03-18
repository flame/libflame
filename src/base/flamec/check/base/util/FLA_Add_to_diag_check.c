
#include "FLAME.h"

FLA_Error FLA_Add_to_diag_check( void* diag_value, FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_null_pointer( diag_value );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

