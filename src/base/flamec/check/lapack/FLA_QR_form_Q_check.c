
#include "FLAME.h"

FLA_Error FLA_QR_form_Q_check( FLA_Obj A, FLA_Obj t )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( t );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

