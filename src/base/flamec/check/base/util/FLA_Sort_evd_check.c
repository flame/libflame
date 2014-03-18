
#include "FLAME.h"

FLA_Error FLA_Sort_evd_check( FLA_Direct direct, FLA_Obj l, FLA_Obj V )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_direct( direct );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( l );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( l );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( V );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( V );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( l, V );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_object_length_equals( V, FLA_Obj_vector_dim( l ) );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

