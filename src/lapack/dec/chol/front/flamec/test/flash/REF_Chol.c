
#include "FLAME.h"

FLA_Error REF_Chol( FLA_Uplo uplo, FLA_Obj A )
{
  return FLA_Chol_blk_external( uplo, A );
}

