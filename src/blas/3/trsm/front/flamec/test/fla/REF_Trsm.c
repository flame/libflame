
#include "FLAME.h"

FLA_Error REF_Trsm( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Trsm_external( side, uplo, trans, diag, alpha, A, B );

  return 0;
}

