
#include "FLAME.h"

FLA_Error FLA_Trsvsx_external( FLA_Uplo uplo, FLA_Trans transa, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y ) 
{
  FLA_Obj x_copy;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsvsx_check( uplo, transa, diag, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, x, &x_copy );
  
  FLA_Copy_external( x, x_copy );

  FLA_Trsv_external( uplo, transa, diag, A, x_copy );

  FLA_Scal_external( beta, y );
  
  FLA_Axpy_external( alpha, x_copy, y );

  FLA_Obj_free( &x_copy );

  return FLA_SUCCESS;
}

