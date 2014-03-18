
#include "FLAME.h"

FLA_Error FLA_Trsmsx( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  return FLA_Trsmsx_external( side, uplo, trans, diag, alpha, A, B, beta, C );
}

