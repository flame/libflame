
#include "FLAME.h"

FLA_Error REF_Syrk( FLA_Uplo uplo, FLA_Trans trans, FLA_Obj alpha, FLA_Obj A, FLA_Obj beta, FLA_Obj C )
{
  FLA_Syrk_external( uplo, trans, alpha, A, beta, C );

  return 0;
}
