
#include "FLAME.h"

FLA_Error FLA_Hevd_check( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj e )
{
  FLA_Error e_val;

  e_val = FLA_Check_valid_evd_type( jobz );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_valid_uplo( uplo );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_floating_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( e );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( A, e );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_square( A );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_vector_dim( e, FLA_Obj_length( A ) );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_col_storage( e );
  FLA_Check_error_code( e_val );
  
  return FLA_SUCCESS;
}

