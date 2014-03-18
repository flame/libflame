
#include "FLAME.h"

FLA_Error FLA_Obj_extract_real_part_check( FLA_Obj a, FLA_Obj b )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( a );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( b );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( b );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( a, b );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( a, FLA_Obj_vector_dim( b ) );
  FLA_Check_error_code( e_val );  

  return FLA_SUCCESS;
}


