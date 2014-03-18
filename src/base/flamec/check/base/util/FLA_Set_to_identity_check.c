
#include "FLAME.h"

FLA_Error FLA_Set_to_identity_check( FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_object_datatype( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

