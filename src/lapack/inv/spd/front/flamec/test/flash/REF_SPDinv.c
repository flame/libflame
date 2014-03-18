
#include "FLAME.h"

FLA_Error REF_SPDinv( FLA_Uplo uplo, FLA_Obj A )
{
  return FLA_SPDinv_blk_external( uplo, A );
}

