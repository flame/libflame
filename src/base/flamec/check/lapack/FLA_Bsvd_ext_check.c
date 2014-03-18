
#include "FLAME.h"

FLA_Error FLA_Bsvd_ext_check( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e,
                              FLA_Obj G, FLA_Obj H,
                              FLA_Svd_type jobu, FLA_Obj U,
                              FLA_Svd_type jobv, FLA_Obj V,
                              FLA_Bool apply_Uh2C, FLA_Obj C )
{
  FLA_Error e_val = FLA_SUCCESS;
  
  FLA_Bsvd_check( uplo, d, e, G, H, jobu, U, jobv, V );

  if ( apply_Uh2C != FALSE ) 
  {
    FLA_Check_identical_object_datatype( U, C );
    FLA_Check_error_code( e_val );

    FLA_Check_object_length_equals( C, FLA_Obj_length( U ) );
    FLA_Check_error_code( e_val );
  }

  return FLA_SUCCESS;
}

