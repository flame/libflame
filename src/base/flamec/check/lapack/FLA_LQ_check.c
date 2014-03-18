
#include "FLAME.h"

FLA_Error FLA_LQ_check( FLA_Obj A, FLA_Obj t )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_vector( t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_storage( t );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( t, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

