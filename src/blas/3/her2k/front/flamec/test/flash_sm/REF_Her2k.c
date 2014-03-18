
#include "FLAME.h"

FLA_Error REF_Her2k( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Her2k_external( uplo, trans, alpha, A, B, beta, C );

  return 0;
}
