
#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_check( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T )
{
  FLA_Error    e_val;

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_storev( storev );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( A, tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_consistent_object_datatype( A, T );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_vector( tau );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( tau, FLA_Obj_min_dim( A ) );
  FLA_Check_error_code( e_val );

  // This is not valid.
  // The width of T can match to either length of width of A, 
  // which depends on how house-holder vectors are accumulated. 
  // e_val = FLA_Check_object_width_equals( T, FLA_Obj_min_dim( A ) );
  // FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

