
#include "FLAME.h"

FLA_Error FLA_Norm1_check( FLA_Obj A, FLA_Obj norm )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, norm );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( norm );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

