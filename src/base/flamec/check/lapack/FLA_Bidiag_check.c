
#include "FLAME.h"

FLA_Error FLA_Bidiag_check( FLA_Obj A, FLA_Obj tu, FLA_Obj tv )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, tu );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_datatype( A, tv );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_vector( tu );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_storage( tu );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_vector( tv );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_storage( tv );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim( tu, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim( tv, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

