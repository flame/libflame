
#include "FLAME.h"

FLA_Error REF_Symm( FLA_Side side, FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Symm_external( side, uplo, alpha, A, B, beta, C );

  return 0;
}
