
#include "FLAME.h"

FLA_Error REF_Trsv( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x )
{
  FLA_Trsv_external( uplo, trans, diag, A, x );

  return 0;
}

