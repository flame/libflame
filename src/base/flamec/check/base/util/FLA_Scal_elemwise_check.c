
#include "FLAME.h"

FLA_Error FLA_Scal_elemwise_check( FLA_Trans trans, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_trans( trans );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, B );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_conformal_dims( trans, A, B );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

