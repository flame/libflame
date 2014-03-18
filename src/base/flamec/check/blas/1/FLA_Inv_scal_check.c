
#include "FLAME.h"

FLA_Error FLA_Inv_scal_check( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  if ( FLA_Obj_is_real( A ) )
  {
    e_val = FLA_Check_consistent_object_datatype( A, alpha );
    FLA_Check_error_code( e_val );
  }
  else
  {
    e_val = FLA_Check_identical_object_precision( A, alpha );
    FLA_Check_error_code( e_val );
  }

  e_val = FLA_Check_if_scalar( alpha );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_divide_by_zero( alpha );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

