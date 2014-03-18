
#include "FLAME.h"

FLA_Error FLA_SPDinv_blk_external( FLA_Uplo uplo, FLA_Obj A )
{
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Error e_val;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_SPDinv_check( uplo, A );

  e_val = FLA_Chol_blk_external( uplo, A );

  if ( e_val != FLA_SUCCESS )
    return e_val;

  e_val = FLA_Trinv_blk_external( uplo, FLA_NONUNIT_DIAG, A );

  if ( e_val != FLA_SUCCESS )
    return e_val;

  FLA_Ttmm_blk_external( uplo, A );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return FLA_SUCCESS;
}

FLA_Error FLA_SPDinv_blk_ext( FLA_Uplo uplo, FLA_Obj A )
{
  return FLA_SPDinv_blk_external( uplo, A );
}
