
#include "FLAME.h"

FLA_Error FLA_LU_piv_check( FLA_Obj A, FLA_Obj p )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_col_vector( p );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_vector_dim_min( p, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

